Readme R-Pakete bauen
---------------------
0. Dokumentation lesen
- "Making R Packages Under Windows.pdf" lesen
- R-Dokumenation "Writing R Extensions”

1. Download the set of “Rtools” or Unix utilities: http://www.murdoch-sutherland.com/Rtools/tools.zip
1. ActiveState Perl installieren

http://www.murdoch-sutherland.com/Rtools/installer.html

2. Paket checken mit

cd C:\Projects\R

R CMD check tmvtnorm
bzw.
Rcmd check tmvtnorm

3. Paket bauen mit 

R CMD build tmvtnorm

4. Paket installieren mit

R CMD INSTALL <package>
R CMD INSTALL tmvtnorm

5. Installierbares .zip Paket aus dem tar.gz bauen
R CMD INSTALL --build tmvtnorm_0.5-1.tar.gz