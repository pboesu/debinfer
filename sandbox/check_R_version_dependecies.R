#hacky script to figure out depends onf imports & depends to correctly set the minimum required R version

#get metadata from CRAN
pdb <- tools::CRAN_package_db()


#get imported and depends packages for debinfer
debinfer_imports <- unlist(strsplit(subset(pdb, Package == 'deBInfer')$Imports, '(,[ \n])'))
debinfer_depends <- unlist(strsplit(subset(pdb, Package == 'deBInfer')$Depends, ', '))

#get depends of those
depends_depends <- subset(pdb, Package %in% debinfer_depends)$Depends
imports_depends <- subset(pdb, Package %in% debinfer_imports, select = c('Package', 'Depends'))

#pull out maximum R requirements
max(unlist(stringr::str_extract_all(depends_depends, '(?<=R[ (>=]{1,6})[0-9.]+')))
max(unlist(stringr::str_extract_all(imports_depends, '(?<=R[ (>=]{1,6})[0-9.]+')))

