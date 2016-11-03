application "fulton";
$C = new Cone(INPUT_RAYS=>[[1,0,0],[1,1,0],[1,2,1],[1,2,2],[1,1,2],[1,0,1]]);
$tv = new NormalToricVariety($C);
$ccoeff = new Vector<Integer>(-1,-1,-1,-1,-1,-1);
$canonical = $tv->add("DIVISOR", COEFFICIENTS=>$ccoeff);
print $canonical->MODULE_GENERATORS; # There is just one, because gorenstein.
$d1coeff = new Vector<Integer>(-3,0,0,0,0,0);
$d1 = $tv->add("DIVISOR", COEFFICIENTS=>$d1coeff);
print $tv->singular_exti_dimension(1,$d1, $d1); # 0
print $tv->singular_exti_dimension(2,$d1, $d1); # 3
$d2coeff = new Vector<Integer>(0,-4,0,0,0,0);
$d2 = $tv->add("DIVISOR", COEFFICIENTS=>$d2coeff);
$cmind2 = $tv->add("DIVISOR", COEFFICIENTS=>$ccoeff - $d2coeff);
$cmind1 = $tv->add("DIVISOR", COEFFICIENTS=>$ccoeff - $d1coeff);
print $tv->singular_exti_dimension(4, $d1, $cmind2);
print $tv->singular_exti_dimension(4, $d2, $cmind1);


print $tv->singular_exti_dimension(1, $d1, $cmind2);
print $tv->singular_exti_dimension(1, $d2, $cmind1);
print $tv->singular_exti_dimension(2, $d1, $cmind2);
print $tv->singular_exti_dimension(2, $d2, $cmind1);
print $tv->singular_exti_dimension(3, $d1, $cmind2);
print $tv->singular_exti_dimension(3, $d2, $cmind1);


