# Use latex to compile to DVI format
$latex = 'latex %O %S';

# Specify the DVIPS step to convert DVI to PS
$dvips = 'dvips %O -o %D %S';

# Specify the PS2PDF step to convert PS to PDF
$ps2pdf = 'ps2pdf %O %S %D';

# Tell latexmk to use the DVI-PS-PDF chain
$pdf_mode = 3;

# Clean up additional files (optional)
$clean_ext = "ps dvi";
