####################################################################################################
## 
## Connecting Bayesian networks and FBA for prediction of metabolic changes
## Methods for estimation of the rates and sampling, as well as obtianing the basal reaction rates
## of control and AD groups are based on the following work:
##
## Gavai AK, Supandi F, Hettling H, Murrell P, Leunissen JA, van Beek JH (2015). 
## Using Bioconductor Package BiGGR for Metabolic Flux Estimation Based on Gene Expression Changes in Brain. 
## PLoS ONE, 10(3), e0119016.
##
####################################################################################################

library("BiGGR")
library("org.Hs.eg.db")
library("BayesNetBP")
library("ggplot2")
load("./data/temp.rda")

set.seed(7654321)
rates.df <- read.csv(file="./data/rates.csv", header = TRUE, row.names = 1)

###################################################
## Set the group: control or AD
###################################################
hf <- 10; pdk <- 0.9;
group <- "control" #  options: "control", "AD"

###################################################
### Build SBML model
###################################################

file.name <- system.file("extdata", "brainmodel_reactions.txt", 
                         package="BiGGR")
reaction.ids <- scan(file.name, what=" ")
reaction.ids <- setdiff(reaction.ids, "R_GLNtN1")
data("H.sapiens_Recon_1")
sbml.model <- buildSBMLFromReactionIDs(reaction.ids, H.sapiens_Recon_1)
rex <- sbml.model@reactions$R_AKGMALtm
rex@reversible <- FALSE
rex@kineticLaw@parameters$LOWER_BOUND@value <- 0
sbml.model@reactions$R_AKGMALtm <- rex

##following term is to be maximized
maximize <- "R_ATPS4m - R_NDPK1m -R_HEX1 - R_PFK - R_PGK + R_PYK"
##specify the external metabolites of the system
externals <- c("M_glc_DASH_D_e", "M_lac_DASH_L_e", "M_ala_DASH_L_e",
               "M_gln_DASH_L_e", "M_h2o_e", "M_co2_e",
               "M_o2_e", "M_h_e", "M_o2s_m",
               "M_adp_c", "M_atp_c", "M_pi_c",
               "M_h_c", "M_nadp_c", "M_nadph_c",
               "M_na1_c", "M_na1_e", "M_gln_DASH_L_c",
               "M_nh4_c", "M_pyr_e")
externals <- setdiff(externals, c("M_gln_DASH_L_e", "M_na1_c", "M_na1_e"))

###################################################
### Get the rates from the initial model
###################################################

##load lying-tunell data
data(lying.tunell.data)

##set equality constraints
equation.vars <- c("R_GLCt1r", "R_L_LACt2r", "R_PYRt2r", "R_GLUDC", "R_G6PDH2r")

equation.values <- c(as.character(
  lying.tunell.data[c("glucose", "lactate", "pyruvate"), "median"]), 
  "R_GLCt1r * 0.32", "R_GLCt1r * 0.069" )

eqns <- list(equation.vars, equation.values)
limfile.path.0 <- tempfile()
createLIMFromSBML(sbml.model, maximize, equations=eqns, 
                  externals=externals, file.name=limfile.path.0)

rates <- getRates(limfile.path.0)
rx.sq <- names(rates)

if (group=="AD") {
  rates <- rates.df$AD ## AD model
} else {
  rates <- rates.df$CTR ## Control model
}
names(rates) <- row.names(rates.df)
rates <- rates[rx.sq]

###################################################
### sample_flux_ensemble
###################################################
uncertain.vars <- data.frame(var=c(equation.vars[c(1,2,3)]),
                             value=as.numeric(c(equation.values[c(1,2,3)])),
                             sd=lying.tunell.data[c("glucose", 
                                                    "lactate", # "glutamine", 
                                                    "pyruvate"), "sd"])
limfile.path.ens <- tempfile()
##Create new LIM model 
equations <- list(c("R_G6PDH2r", "R_GLUDC", "R_G3PD2m") , 
                  c("R_GLCt1r * 0.069", "R_GLCt1r * 0.32", "0"))
createLIMFromSBML(sbml.model, maximize, externals=externals, 
                  file.name=limfile.path.ens, equations=equations) 

##sample feasible flux distributions with MCMC
ensemble <- sampleFluxEnsemble(limfile.path.ens, uncertain.vars, 
                               x0=rates, iter=1e5, burninlength=1e4, 
                               outputlength=1e4, type="mirror", jmp=0.1)

atp.prod.ens <- eval(parse(text=maximize), envir=data.frame(ensemble))
mean(atp.prod.ens)

###################################################
### compute fold-change
###################################################

###### select control or AD Bayesian network
if (group=="AD") {
  tree.init.p <- tree.init.ad ## AD model
} else {
  tree.init.p <- tree.init.ctrl ## Control model
}
enzymes <- unique(unlist(constraint.rx)) # enzymes in reations

# perturb model by changing HIF1A levels
tree.post <- AbsorbEvidence(tree.init.p, c("HIF1A"), list(hf))
## compute mean values of HIF1A
SummaryMarginals(Marginals(tree.init.p, "HIF1A"))$Mean

exp.0 <- SummaryMarginals(Marginals(tree.init.p, enzymes))$Mean # mean expression from initial model
exp.1 <- SummaryMarginals(Marginals(tree.post, enzymes))$Mean # mean expression from perturbed model
names(exp.0) <- enzymes
names(exp.1) <- enzymes

fold.gene <- exp.1/exp.0 # gene expression fold change
# compute reaction fold change by taking average
fold.rx <- c() 
for(i in 1:length(constraint.rx)) {
  fold.rx[i] <- mean(fold.gene[constraint.rx[[i]]])
}
rxs <- names(constraint.rx)
names(fold.rx) <- rxs
rxs[10] <- "R_SSALxm"

## compute the additional constraints by rate*fold-change
constraints <- rates[rxs] * fold.rx
constraints["R_PDHm"] <- constraints["R_PDHm"] * pdk  #########################################

## variances of flux
esbl.sub <- ensemble[, rxs]
var.flux <- apply(esbl.sub, 2, var)

###################################################
### Include additional constraints into model
###################################################

ind <- 1:10 # there are 10 reactions as new constraints, here we add them all
equation.vars.2 <- c(equation.vars[3:5], rxs[ind])

if (group=="AD") {
  ## AD model, the rates of externals are updated by the AD rates accordingly
  equation.values.2 <- c("-0.0439", "R_GLCt1r * 0.8978", "R_GLCt1r * 0.2099", as.character(constraints[ind])) 
} else {
  ## Control model
  equation.values.2 <- c(equation.values[3:5], as.character(constraints[ind])) 
}

eqns.2 <- list(equation.vars.2, equation.values.2)
limfile.path <- tempfile()
createLIMFromSBML(sbml.model, maximize, equations=eqns.2, 
                  externals=externals, file.name=limfile.path)
lim.model <- Setup(limfile.path) # extract matrices
rates.3 <- lsei(A=lim.model$A[76:85,], B=lim.model$B[76:85], 
                E=lim.model$A[c(1:75),], F=lim.model$B[c(1:75)], 
                G=lim.model$G, H=lim.model$H)$X
names(rates.3) <- names(rates)

###################################################
### Sampling
###################################################

lim.model$VarB <- var.flux
ensbl <- xsample(A=lim.model$A[76:85,], B=lim.model$B[76:85], 
                 E=lim.model$A[c(1:75),], F=lim.model$B[c(1:75)], 
                 G = lim.model$G, H = lim.model$H, sdB = sqrt(lim.model$VarB), W = 1, 
                 iter = 10000, outputlength = 10000, burninlength = 2000, 
                 type = "mirror", jmp = 5, tol = sqrt(.Machine$double.eps), 
                 x0 = NULL, fulloutput = FALSE, test = TRUE)

sample.x <- data.frame(ensbl$X)
names(sample.x) <- lim.model$Unknowns

rates.ens <- colMeans(sample.x)
rates.var <- apply(sample.x, 2, var)
rates.sd <- sqrt(rates.var)

rates.output <- data.frame(rates.ens, rates.sd)
f.name <- paste0(group, "_", as.character(hf*10), ".csv")
# write.csv(rates.output, file=f.name)
rxnm <- rownames(rates.output)

## compute ATP production
atp <- eval(parse(text=maximize), envir=data.frame(sample.x))
mean(atp) 
## mean ATP level @ HIF 8~9.5~13: 
## control 6.64~8.05~11.38; 
## AD 6.09~6.42~7.23;  ~4.27~4.48
## model without HIF1A perturbation 9.21.

#######################################################################
## plot pathway
## Glycolysis
relevant.species <- c("M_glc_DASH_D_c", "M_g6p_c", "M_f6p_c",
                      "M_fdp_c", "M_dhap_c", "M_g3p_c",
                      "M_13dpg_c", "M_3pg_c", "M_2pg_c",
                      "M_pep_c", "M_pyr_c",
                      "M_6pgl_c", "M_6pgc_c", "M_ru5p_DASH_D_c", 
                      "M_xu5p_DASH_D_c", "M_r5p_c", "M_g3p_c", "M_s7p_c",
                      "M_lac_DASH_L_e", "M_lac_DASH_L_c", "M_pyr_m", "M_pyr_c")

relevant.reactions <- c("R_HEX1", "R_PGI", "R_PFK", "R_FBA", "R_TPI", 
                        "R_GAPD", "R_PGK", "R_PGM", "R_ENO", "R_PYK",
                        "R_G6PDH2r", "R_PGL", "R_GND", "R_RPE", "R_RPI", "R_TKT1",
                        "R_L_LACt2r", "R_LDH_L", "R_PDHm", "R_PYRt2m")

hd.ens <- sbml2hyperdraw(sbml.model, rates=rates.ens, 
                         relevant.species=relevant.species, 
                         relevant.reactions=relevant.reactions,
                         layoutType="dot", plt.margins=c(20, 0, 20, 80))
plot(hd.ens)

## TCA cycle
relevant.species <- c("M_cit_m", "M_icit_m" , "M_akg_m",
                      "M_succoa_m", "M_succ_m", "M_fum_m",
                      "M_mal_DASH_L_m", "M_oaa_m")
relevant.reactions <- c("R_CSm", "R_ACONTm", "R_ICDHxm",
                        "R_AKGDm", "R_SUCOAS1m", "R_SUCD1m",
                        "R_FUMm", "R_MDHm", "R_ICDHyrm", "R_ME1m",
                        "R_ME2m", "R_ASPTAm","R_AKGMALtm", "R_GLUDym",
                        "R_ABTArm", "R_SSALxm","R_CITtam")
hd.ens.2 <- sbml2hyperdraw(sbml.model, rates=rates.ens,
                           relevant.reactions=relevant.reactions,
                           relevant.species=relevant.species,
                           layoutType="circo", plt.margins=c(150, 235, 150, 230))

plot(hd.ens.2)
#######################################################################