#!/usr/bin/env perl

$infile = $ARGV[0];

chomp($xsize = `mincinfo -dimlength xspace $infile`);
chomp($ysize = `mincinfo -dimlength yspace $infile`);

print STDOUT "Sizes: $xsize:$ysize\n";

$tmp = `mincextract -float -unsigned $infile`;

$code = 'f' . $xsize*$ysize;
@data = unpack($code, $tmp);

$min = $max = $data[0];
foreach (@data){
   if($_ > $max){
      $max = $_;
      }
   
   if($_ < $min){
      $min = $_;
      }
   }
print STDOUT "Min: $min   Max: $max\n";

# play trajectories
$num = 60000;

# simple "ray"
#for($c=0; $c<$num; $c++){
#   $value = $data[($c*$xsize)+$c];
#   print "[$c] $value\n";
#   push(@out, $value);
#   push(@coords, "$c,$c,0");
#   }

# random points
for($c=0; $c<$num; $c++){
   
   $x = int(rand($xsize-1));
   $y = int(rand($ysize-1));
   
   $value = $data[($y*$xsize)+$x];
#   print "[$c] [$x:$y] $value\n";
   push(@out, $value);
   push(@coords, "$x,$y,0");
   }

# write out data stage
open(FH, ">raw_data");
$tmp = pack("f$num", @out);
syswrite(FH, $tmp);
close(FH);

open(FH, ">coord_data");
foreach (@coords){
   print FH "$_\n";
   }
close(FH);
