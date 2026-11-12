@ARGV = grep { !/^-{1,2}(outdir|output-directory|auxdir|aux-directory)=/ } @ARGV;
$pdf_mode = 1;
$aux_dir = $out_dir = '.build';
$success_cmd = '[ "%W" != "." ] && cp -f %W/%R.pdf ./';
