#load('src/linnaeus-scripts/collapsibleTree/sand/Z2_Ltree_pie.Robj')
# load('~/Dropbox/PJ/linnaeus_manuscript/collapsibleTrees/Trees/A7_Ltree_pie.Robj')
load('~/Dropbox/scartrace manuscript/collapsibleTrees/Trees/A7_Ltree_pie.Robj')

f_zoom2 = !is.na(linnaeus.colors_adult$zoom2)
ctypes_exo = paste(linnaeus.colors_adult$Cell.type[f_zoom2], "pancreas", "exo")
ctypes_endo = paste(linnaeus.colors_adult$Cell.type[f_zoom2], "pancreas", "endo")
col_pan = c(colorRampPalette(c("#0029b2","#b4dbfd"))(9), "#007a00",  colorRampPalette(c("#b20500","#feebed"))(9), "#930093")

ct_pan = c(ctypes_endo, ctypes_exo)

ctypes = ct_pan
#ctypes =  grep("alpha|beta|delta|duct", x=ctypes, perl=T, value = T )
ct_colors = col_pan

orit = Clone(LINNAEUS.pie)

t <- data.tree::Traverse(orit, "level")
sapply(t, function(x){if(grepl("endo", x$name)){x$Cell.type = paste(x$Cell.type, "endo")}})
sapply(t, function(x){if(grepl("exo", x$name)){x$Cell.type = paste(x$Cell.type, "exo")}})

get_pieNode(orit, ctypes=ctypes)
#ttt = Clone(orit$nd0_27$nd0_27_1)
ttt = Clone(orit)

namess = ttt$Get(function(x) if(x$isScar) x$name)
namess = namess[!is.na(namess)]
names(namess) = namess

