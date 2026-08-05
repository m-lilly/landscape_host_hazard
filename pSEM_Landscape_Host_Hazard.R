
#######Piecewise Structural Equation Modeling code for manuscript titled, 'Urbanization shifts roles of mammalian hosts in tick-borne hazard emergence'########

#Read data
data <-read.csv("pSEM_data.csv", header = T)

#Load libraries
library(piecewiseSEM)
library(glmmTMB)
library(lme4)

############## SEM PIECES #################

###Final glmm model pieces####
m_iscap <- glmmTMB(ISCAP_Larva ~ Sex_01 + DON_IS + (1|Site), family = nbinom2(link="log"),
                    data = data)

m_don <- glmmTMB(DON_IS ~  Connect_site_1000 + (1|Park), family = gaussian(),
                 data = data)

m_din <- glmmTMB(DIN_Borrelia ~ DON_IS + Connect_site_1000*SM_abundance +(1|Site), family = gaussian(),
                 data = data)

m_bb <- glm(Borrelia_qPCR ~ DIN_Borrelia + Weight + Connect_site_1000*MM_abundance, family = binomial(link="logit"),
            data = data)
summary(m_bb)

m_bab <- glmer(Babesia_qPCR ~ Borrelia_qPCR  + Weight + DON_IS + (1|Park), offset = (log(Effort_qPCR)), family = binomial(link="logit"),
               data = data)

mm_ra <- glmmTMB(MM_abundance ~ Impervious_site + (1|Site), family= gaussian(),
                 data = data)

sm_ra <- glmmTMB(SM_abundance ~ Connect_site_1000 + (1|Site), family= gaussian(),
                 data = data)

##Test for spatial autocorrelation
library(spdep)
sm.coords <- cbind(data$Longitude, SM_subset$Latitude)
sm.5nn <- knearneigh(sm.coords, k=5, longlat = TRUE)
sm.5nn.nb <- knn2nb(sm.5nn)
plot(sm.5nn.nb, sm.coords)
listw <- nb2listw(sm.5nn.nb)

#Moran's I
moran.mc(residuals(m_iscap), listw, 1000, zero.policy = T)
moran.mc(residuals(m_don), listw, 1000, zero.policy = T)
moran.mc(residuals(m_din), listw, 1000, zero.policy = T)
moran.mc(residuals(m_bb), listw, 1000, zero.policy = T)
moran.mc(residuals(m_bab), listw, 1000, zero.policy = T)
moran.mc(residuals(mm_ra), listw, 1000, zero.policy = T)
moran.mc(residuals(sm_ra), listw, 1000, zero.policy = T)


###piecewise SEM model####
model_full <- psem(m_don, m_bb, m_din, m_bab, mm_ra, sm_ra,m_iscap,
                   (SM_abundance %~~% DON_IS),  (SM_abundance %~~% MM_abundance), 
                   (MM_abundance %~~% DIN_Borrelia), (Babesia_qPCR%~~% Connect_site_1000),
                   (MM_abundance%~~% Connect_site_1000), (MM_abundance%~~% DON_IS), (ISCAP_Larva%~~% MM_abundance),
                   (ISCAP_Larva%~~% DIN_Borrelia), (DON_IS%~~% Sex_01),  (DON_IS%~~%Impervious_site))
summary(model_full, standardize = "scale", conserve = TRUE)



#####Calculating indirect effects to Borrelia and Babesia mouse infection status ######
# --- Coefficients and standard errors from summary(model_full) ---
coef_table <- summary(model_full)$coefficients
rownames(coef_table) <- paste(coef_table$Response, coef_table$Predictor, sep="~")

# ----------------------------
# Path 1: Connect_site_1000 → DON_IS → DIN_Borrelia → Borrelia_qPCR
# ----------------------------

# Step 1: Connect_site_1000 → DON_IS
a <- coef_table["DON_IS~Connect_site_1000", "Estimate"]
SE_a <- as.numeric(coef_table["DON_IS~Connect_site_1000", "Std.Error"])

# Step 2: DON_IS → DIN_Borrelia
b <- coef_table["DIN_Borrelia~DON_IS", "Estimate"]
SE_b <- as.numeric(coef_table["DIN_Borrelia~DON_IS", "Std.Error"])

# Step 3: DIN_Borrelia → Borrelia_qPCR
c <- coef_table["Borrelia_qPCR~DIN_Borrelia", "Estimate"]
SE_c <- as.numeric(coef_table["Borrelia_qPCR~DIN_Borrelia", "Std.Error"])


# --- 3. Compute indirect effect ---
IE <- a * b * c

#IE <- a * b

# --- 4. Compute standard error for three-step path ---
SE_IE <- sqrt((b*c)^2 * (SE_a)^2 + (a*c)^2 * (SE_b)^2 + (a*b)^2 * (SE_c)^2 )

#SE_IE <- sqrt( (b^2 * SE_a^2) + (a^2 * SE_b^2) )

# --- 5. Compute z-value and p-value ---
z_val <- IE / SE_IE
p_val <- 2 * (1 - pnorm(abs(z_val)))

# --- 6. Print results ---
IE
SE_IE
z_val
p_val

##p_val = 0.07 so not quite a significant effect of connectivity on Borrelia qPCR
#Yes significant indirect effect on DIN... but also significant direct effect


# ----------------------------
# Path 2: Impervious_site → MM_abundance → Borrelia_qPCR
# ----------------------------
a2 <- coef_table["MM_abundance~Impervious_site", "Estimate"]
SE_a2 <- as.numeric(coef_table["MM_abundance~Impervious_site", "Std.Error"])

b2 <- coef_table["Borrelia_qPCR~MM_abundance", "Estimate"]
SE_b2 <- as.numeric(coef_table["Borrelia_qPCR~MM_abundance", "Std.Error"])

IE2 <- a2 * b2
SE_IE2 <- sqrt( (b2^2 * SE_a2^2) + (a2^2 * SE_b2^2) )
z2 <- IE2 / SE_IE2
p2 <- 2 * (1 - pnorm(abs(z2)))

IE2; SE_IE2; z2; p2

##Not significant: p = 0.08

# ----------------------------
# Path 3: Connect_site_1000 → SM_abundance → DIN_Borrelia → Borrelia_qPCR
# ----------------------------
a3 <- coef_table["SM_abundance~Connect_site_1000", "Estimate"]
SE_a3 <- as.numeric(coef_table["SM_abundance~Connect_site_1000", "Std.Error"])

b3 <- coef_table["DIN_Borrelia~SM_abundance", "Estimate"]
SE_b3 <-  as.numeric(coef_table["DIN_Borrelia~SM_abundance", "Std.Error"])

c3 <- coef_table["Borrelia_qPCR~DIN_Borrelia", "Estimate"]
SE_c3 <-  as.numeric(coef_table["Borrelia_qPCR~DIN_Borrelia", "Std.Error"])

IE3 <- a3 * b3 * c3
SE_IE3 <- sqrt( (b3*c3)^2 * SE_a3^2 + (a3*c3)^2 * SE_b3^2 + (a3*b3)^2 * SE_c3^2 )
z3 <- IE3 / SE_IE3
p3 <- 2 * (1 - pnorm(abs(z3)))

IE3; SE_IE3; z3; p3
#p = 0.09

# ----------------------------
# Path 3.5: Connect_site_1000 → SM_abundance → DIN_Borrelia → Borrelia_qPCR
# ----------------------------

IE3 <- b3 * c3
SE_IE3 <- sqrt( b3^2 * SE_c3^2 + c3^2 * SE_b3^2)
z3 <- IE3 / SE_IE3
p3 <- 2 * (1 - pnorm(abs(z3)))

IE3; SE_IE3; z3; p3
#p = also not significant indirect from small RA to Borrelia 
#
# ----------------------------
# Path 4: Connect_site_1000 → DON_IS → ISCAP_Larva
# ----------------------------
a4 <- coef_table["DON_IS~Connect_site_1000", "Estimate"]
SE_a4 <- as.numeric(coef_table["DON_IS~Connect_site_1000", "Std.Error"])

b4 <- coef_table["ISCAP_Larva~DON_IS", "Estimate"]
SE_b4 <-  as.numeric(coef_table["ISCAP_Larva~DON_IS", "Std.Error"])

IE4 <- a4 * b4
SE_IE4 <- sqrt( (b4)^2 * SE_a4^2 + (a4)^2 * SE_b4^2)
z4 <- IE4 / SE_IE4
p4 <- 2 * (1 - pnorm(abs(z4)))

IE4; SE_IE4; z4; p4
#p = 0.02, IE = 0.16
#####
