application "fulton";

$c = new Cone(INPUT_RAYS=>[[1,0,1],[1,1,0],[1,2,0],[1,1,1],[1,0,2]]);
$dc = new Cone(INPUT_RAYS=>$c->FACETS);
$tv = new NormalToricVariety($dc);
print $tv->ISOLATED; # Is isolated
print $tv->GORENSTEIN; # Is Gorenstein
print $tv->WEIGHT_CONE->HILBERT_BASIS; # Has very few gens for algebra
$canonicalCoeff = new Vector(-1,-1,-1,-1);
$div1Coeff = new Vector(7,-1,-1,-1);
$div1 = $tv->add("DIVISOR", COEFFICIENTS=>$div1Coeff);
$canonical = $tv->add("DIVISOR", COEFFICIENTS=>$canonicalCoeff);
print $tv->singular_exti_dimension(1,$div1, $div1);
print $tv->singular_exti_dimension(2,$div1, $div1);
print $tv->singular_exti_dimension(3,$div1, $div1);
print $tv->singular_exti_dimension(4,$div1, $div1);
$div2Coeff = new Vector(-5,0,0,0);
$div2 = $tv->add("DIVISOR", COEFFICIENTS=>$div2Coeff);
$kmd1 = $tv->add("DIVISOR", COEFFICIENTS=>$canonicalCoeff - $div1Coeff);
$kmd2 = $tv->add("DIVISOR", COEFFICIENTS=>$canonicalCoeff - $div2Coeff);
$d1mk = $tv->add("DIVISOR", COEFFICIENTS=> - $canonicalCoeff + $div1Coeff);
$d2mk = $tv->add("DIVISOR", COEFFICIENTS=> - $canonicalCoeff + $div2Coeff);
print $tv->singular_exti_dimension(2, $div1, $kmd2);
print $tv->singular_exti_dimension(2, $div2, $kmd1);

print $tv->singular_exti_dimension(3, $div2, $kmd1);
print $tv->singular_exti_dimension(3, $div1, $kmd2);

print $tv->singular_exti_dimension(4, $div1, $kmd2);
print $tv->singular_exti_dimension(4, $div1, $d2mk);
print $tv->singular_exti_dimension(4, $div2, $kmd1);
print $tv->singular_exti_dimension(4, $div2, $d1mk);
print $tv->singular_tori_dimension(1, $div2, $div1);

print $tv->singular_exti_dimension(5, $div1, $kmd2);
print $tv->singular_exti_dimension(5, $div1, $d2mk);
print $tv->singular_exti_dimension(5, $div2, $kmd1);
print $tv->singular_exti_dimension(5, $div2, $d1mk);
print $tv->singular_tori_dimension(2, $div2, $div1);

print $tv->singular_exti_dimension(5, $div1, $kmd2);
print $tv->singular_exti_dimension(5, $div2, $kmd1);
print $tv->singular_tori_dimension(2, $div2, $div1);

print $tv->singular_exti_dimension(6, $div1, $kmd2);
print $tv->singular_exti_dimension(6, $div2, $kmd1);
print $tv->singular_tori_dimension(3, $div2, $div1);

print $tv->singular_exti_dimension(7, $div1, $kmd2);
print $tv->singular_exti_dimension(7, $div2, $kmd1);
print $tv->singular_tori_dimension(4, $div2, $div1);

print $tv->singular_exti_dimension(8, $div1, $kmd2);
print $tv->singular_exti_dimension(8, $div2, $kmd1);
print $tv->singular_tori_dimension(5, $div2, $div1);


application "fulton";

$c = new Cone(INPUT_RAYS=>[[1,0,0],[1,0,1],[1,1,0],[1,1,1]]);
$dc = new Cone(INPUT_RAYS=>$c->FACETS);
$tv = new NormalToricVariety($dc);
$canonicalCoeff = new Vector(-1,-1,-1,-1);
$n = 8;
$e = 15;
for(my $i=-$n; $i<=$n; $i++){
   my $div1Coeff = new Vector($i,0,0,0);
   my $div1 = $tv->add("DIVISOR", COEFFICIENTS=>$div1Coeff);
   my $kmd1 = $tv->add("DIVISOR", COEFFICIENTS=>$canonicalCoeff - $div1Coeff);
   for(my $j=-$n; $j<=$n; $j++){
      my $div2Coeff = new Vector($j,0,0,0);
      my $div2 = $tv->add("DIVISOR", COEFFICIENTS=>$div1Coeff);
      my $kmd2 = $tv->add("DIVISOR", COEFFICIENTS=>$canonicalCoeff - $div1Coeff);
      print "i: ",$i," j: ",$j,"\n";
      for(my $k=4; $k<$e; $k++){
         print "Ext^",$k," ";
         my $e12 = $tv->singular_exti_dimension($k, $div1, $kmd2);
         print "12 ";
         my $e21 = $tv->singular_exti_dimension($k, $div1, $kmd2);
         print "21 ";
         my $t = $tv->singular_tori_dimension($k-3, $div2, $div1);
         print "tor\n";
         my $test = ($e12 == $e21) && ($e12 == $t);
         if(!$test){
            die "Failed for $i and $j";
         }
      }
   }
}
