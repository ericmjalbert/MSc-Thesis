set terminal postscript eps color enhanc
set output "converg_total.eps"

set xl "t"
set yl "Average M"

set key outside

plot \
"<paste ./256_total.dat ./128_total.dat" us 1:(abs($2 - $5)/$5) wi lin t "M^{256} - M^{128}",\
"<paste ./512_total.dat ./256_total.dat" us 1:(abs($2 - $5)/$5) wi lin t "M^{512} - M^{256}",\
"<paste ./1024_total.dat ./512_total.dat" us 1:(abs($2 - $5)/$5) wi lin  t "M^{1024} - M^{512}",\
"<paste ./2048_total.dat ./1024_total.dat" us 1:(abs($2 - $5)/$5) wi lin t "M^{2048} - M^{1024}",\
"<paste ./4096_total.dat ./2048_total.dat" us 1:(abs($2 - $5)/$5) wi lin t "M^{4096} - M^{2048}",\
"<paste ./8192_total.dat ./4096_total.dat" us 1:(abs($2 - $5)/$5) wi lin t "M^{8192} - M^{4096}",\
"<paste ./16384_total.dat ./8192_total.dat" us 1:(abs($2 - $5)/$5) wi lin t "M^{16384} - M^{8192}"
