######
# 06 #
####################################################################################################
####################################################################################################
############################################# KeyGeneExp ###########################################
####################################################################################################
####################################################################################################

diffGeneExp <- read.csv("diffGeneExp.txt", sep = "\t")
load("GS_MM.RData")

colnames(diffGeneExp) <- gsub("\\.", "-", colnames(diffGeneExp))



a_greenyellow <- GS_MM[which(GS_MM$moduleColor == "greenyellow"),
                       which((colnames(GS_MM) %in% c("probes", "moduleColor", "GS.", "p.GS.", "MMgreenyellow", "p.MMlightyellow") == TRUE))]
max(a_greenyellow$GS.)
max(a_greenyellow$MMgreenyellow)

intersect(which(a_greenyellow$GS. >= 0.3), which(a_greenyellow$MMgreenyellow >= 0.5))

a_tan <- GS_MM[which(GS_MM$moduleColor == "tan"),
               which((colnames(GS_MM) %in% c("probes", "moduleColor", "GS.", "p.GS.", "MMtan", "p.MMtan") == TRUE))]

max(a_tan$GS.)
min(a_tan$GS.)
max(a_tan$MMtan)
min(a_tan$MMtan)
intersect(which(a_tan$GS. <= -0.3), which(a_tan$MMtan >= 0.5))



gene_tan <- GS_MM[intersect(which(a_tan$GS. <= -0.3), which(a_tan$MMtan >= 0.5)), 
                  1]

gene_greenyellow <- GS_MM[intersect(which(a_greenyellow$GS. >= 0.3), which(a_greenyellow$MMgreenyellow >= 0.5)), 
                          1]

keyGeneExp_tan <- diffGeneExp[which((diffGeneExp$ID %in% gene_tan) == TRUE), ]
keyGeneExp_greenyellow <- diffGeneExp[which((diffGeneExp$ID %in% gene_greenyellow) == TRUE), ]

save(keyGeneExp_tan, file = "keyGeneExp_tan.RData")
save(keyGeneExp_greenyellow, file = "keyGeneExp_greenyellow.RData")