PACKAGE_NAME=tmvtnorm

check:
	R CMD check $(PACKAGE_NAME)

build:
	R CMD build $(PACKAGE_NAME)
	R CMD INSTALL --build $(PACKAGE_NAME)_*.tar.gz

install:
	R CMD INSTALL $(PACKAGE_NAME)