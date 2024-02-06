.onAttach <- function(lib, pkg) {     
    
    requireNamespace("fabia", quietly = TRUE)

    version <- packageDescription("miRSM", fields = "Version")
    packageStartupMessage( "Citation: Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD,","\n",
      "miRSM: an R package to infer and analyze miRNA sponge modules in heterogeneous data,","\n",
      "RNA Biol. 2021 Dec;18(12):2308-2320.","\n",
      "BibTex: enter 'toBibtex(citation(\"miRSM\"))'","\n\n",
      "Homepage: https://github.com/zhangjunpeng411/miRSM","\n\n",
      "miRSM Package Version ", version, "\n")
}
