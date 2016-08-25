application "fulton";

$c = new Cone(INPUT_RAYS=>[[1,0,0],[1,0,1],[1,1,0],[2,1,2]]);
$tv = new NormalToricVariety($c);
print $tv->ISOLATED; # Is isolated
print $tv->GORENSTEIN; # Is not Gorenstein
print $tv->WEIGHT_CONE->HILBERT_BASIS; # Has very few gens for algebra
$canonicalCoeff = new Vector(-1,-1,-1,-1);
$div1Coeff = new Vector(-1,0,0,0);
$div1 = $tv->add("DIVISOR", COEFFICIENTS=>$div1Coeff);
print $tv->singular_exti_dimension(1,$div1, $div1);
print $tv->singular_exti_dimension(2,$div1, $div1);
print $tv->singular_exti_dimension(3,$div1, $div1);
print $tv->singular_exti_dimension(4,$div1, $div1);
$div2Coeff = new Vector(0,-1,0,0);
$div2 = $tv->add("DIVISOR", COEFFICIENTS=>$div2Coeff);
$kmd1 = $tv->add("DIVISOR", COEFFICIENTS=>$canonicalCoeff - $div1Coeff);
$kmd2 = $tv->add("DIVISOR", COEFFICIENTS=>$canonicalCoeff - $div2Coeff);
print $tv->singular_exti_dimension(2, $div1, $kmd2);
print $tv->singular_exti_dimension(2, $div2, $kmd1);

print $tv->singular_exti_dimension(3, $div2, $kmd1);
print $tv->singular_exti_dimension(3, $div1, $kmd2);

print $tv->singular_exti_dimension(4, $div1, $kmd2);
print $tv->singular_exti_dimension(4, $div2, $kmd1);
print $tv->singular_tori_dimension(1, $div2, $div1);

print $tv->singular_exti_dimension(5, $div1, $kmd2);
print $tv->singular_exti_dimension(5, $div2, $kmd1);
print $tv->singular_tori_dimension(2, $div2, $div1);

print $tv->singular_exti_dimension(6, $div1, $kmd2);
print $tv->singular_exti_dimension(6, $div2, $kmd1);
print $tv->singular_tori_dimension(3, $div2, $div1);

print $tv->singular_exti_dimension(7, $div1, $kmd2);
print $tv->singular_exti_dimension(7, $div2, $kmd1);
print $tv->singular_tori_dimension(4, $div2, $div1);
