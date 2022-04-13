
load('~/Dropbox/PJ/linnaeus_manuscript/collapsibleTrees/Trees/Z2_Ltree_pie.Robj')
#load('~/Dropbox/PJ/linnaeus_manuscript/collapsibleTrees/Trees/Z4_Ltree_pie.Robj')
#load('~/Dropbox/PJ/linnaeus_manuscript/collapsibleTrees/Trees/Z5_Ltree_pie.Robj')


f_zoom3 = !is.na(linnaeus.larva$zoom3)
ct_colors = linnaeus.larva$color3[f_zoom3]
ctypes = linnaeus.larva$Cell.type[f_zoom3]


orit = Clone(LINNAEUS.pie)

get_pieNode(orit, ctypes=ctypes)
ttt = Clone(orit)



namess = ttt$Get(function(x) if(x$isScar) x$name)
namess = namess[!is.na(namess)]
names(namess) = namess

