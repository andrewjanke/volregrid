#! /usr/bin/env perl


use warnings "all";
use Data::Dumper;

$large = 1.0e100;

# grunge in the info
foreach (@ARGV){
   
   @args = ('mincinfo',
            '-dimlength', 'xspace',
            '-dimlength', 'yspace',
            '-dimlength', 'zspace',
            $_);
   @sizes = split(/\n/, `@args`);
   
#   @args = ('mincinfo',
#            '-attvalue', 'xspace:start',
#            '-attvalue', 'yspace:start',
#            '-attvalue', 'zspace:start',
#            $_);
#   @starts = split(/\n/, `@args`);
#   
#   @args = ('mincinfo',
#            '-attvalue', 'xspace:step',
#            '-attvalue', 'yspace:step',
#            '-attvalue', 'zspace:step',
#            $_);
#   @steps = split(/\n/, `@args`);
#   
#   @args = ('mincinfo',
#            '-attvalue', 'xspace:direction_cosines',
#            '-attvalue', 'yspace:direction_cosines',
#            '-attvalue', 'zspace:direction_cosines',
#            $_);
#   @dircos = split(/\n/, `@args`);
#   
#   ($mat[0][0], $mat[1][0], $mat[2][0]) = split(/\ /, $dircos[0], 3);
#   ($mat[0][1], $mat[1][1], $mat[2][1]) = split(/\ /, $dircos[1], 3);
#   ($mat[0][2], $mat[1][2], $mat[2][2]) = split(/\ /, $dircos[2], 3);
#   
#   # generate each of the eight corners
#   for($i=0; $i<3; $i++){
#      @{$coords[$i]} = ($starts[$i], $starts[$i] + ($steps[$i] * ($sizes[$i]-1)));
#      }
#   
#   # find the min and max's
#   @min = ($large, $large, $large);
#   @max = (-$large, -$large, -$large);
#   foreach $x (@{$coords[0]}){
#      foreach $y (@{$coords[1]}){
#         foreach $z (@{$coords[2]}){
#           
#            @coord = ($x, $y, $z);
#            @res = &v_mult(*coord, *mat);
#            
#            for($i=0; $i<3; $i++){
#               if($res[$i] > $max[$i]){
#                  $max[$i] = $res[$i];
#                  }
#               if($res[$i] < $min[$i]){
#                  $min[$i] = $res[$i];
#                  }
#               }
#            print "BEFORE: @coord    RESULT: @res\n";
#            }
#         }
#      }
   
   # find the min and max's
   @min = ($large, $large, $large);
   @max = (-$large, -$large, -$large);
   foreach $x (0, $sizes[0]-1){
      foreach $y (0, $sizes[1]-1){
         foreach $z (0, $sizes[2]-1){
            
            chomp($tmp = `voxeltoworld $_ $x $y $z`);
            (@res) = split(/\ /, $tmp);          
            for($i=0; $i<3; $i++){
               if($res[$i] > $max[$i]){
                  $max[$i] = $res[$i];
                  }
               if($res[$i] < $min[$i]){
                  $min[$i] = $res[$i];
                  }
               }
            print "VOXEL: $x $y $z   WORLD: @res\n";
            }
         }
      }
   }

print "MIN: @min\n";
print "MAX: @max\n";


$sampling = 0.231729*2;

@args = ('-xstart', $min[0],
         '-ystart', $min[1],
         '-zstart', $min[2],
         '-xstep', $sampling,
         '-ystep', $sampling,
         '-zstep', $sampling,
         '-xlength', int(($max[0] - $min[0]) / $sampling) + 1,
         '-ylength', int(($max[1] - $min[1]) / $sampling) + 1,
         '-zlength', int(($max[2] - $min[2]) / $sampling) + 1,
         '-regrid_radius', $sampling,
         );

print STDOUT "args: " . join(" ", @args) . "\n";
         

# matrix multiplication
sub v_mult {
   local (*vect, *mat) = @_;
#   print Dumper(@mat);
   
   for($i=0; $i<3; $i++){
      $res[$i] = 0;
      for($j=0; $j<3; $j++){
         $res[$i] += $mat[$i][$j] * $vect[$j];
         }
      }
   
   return (@res);
   }
