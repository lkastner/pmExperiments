use Benchmark qw(:all) ;
application "fulton";
prefer "libnormaliz";
$n = 283200000;
$q = 15587947;
print gcd($n,$q);
$cone = new Cone(INPUT_RAYS=>[[0,1],[$n,$q]]);
$cone = new Cone(RAYS=>$cone->RAYS, FACETS=>$cone->FACETS);
$t0 = Benchmark->new; $cone->HILBERT_BASIS; $t1 = Benchmark->new;
print timestr(timediff($t1, $t0));
$cqs = new CyclicQuotient(N=>$n, Q=>$q);
$s0 = Benchmark->new; @test = cqs_to_hb($cqs); $s1 = Benchmark->new;
print timestr(timediff($s1, $s0));
print $cqs->DUAL_CONTINUED_FRACTION->dim;
print $cone->HILBERT_BASIS->rows;
print scalar(@test);

$n = 28320;
$q = 1559;

# Good 3dim examples:
$c = new Cone(INPUT_RAYS=>[[1,0,0],[32,34,51],[3,-29,77]]);
$c = new Cone(INPUT_RAYS=>[[1,0,0],[32,34,51],[3,-29,7]]);
$c = new Cone(INPUT_RAYS=>[[1,0,0],[32,34,51],[3,-9,7]]);

$hb = $c->HILBERT_BASIS;
$hb = new Matrix<Rational>($hb);
$test = new Matrix(threedim_case($c));
foreach my $h (@$hb){
   my $n = grep($_ == $h, @$test);
   print $h,": ",$n,"\n";
}

# C++ code
$dcf = new Vector<Integer>([2]);
$hb = new Matrix<Integer>([[0,1],[1,1],[2,1]]);
print module_generator_recursion(1,1,$dcf,$hb);
print module_generator_recursion(2,1,$dcf,$hb);
print module_generator_recursion(1,2,$dcf,$hb);
print module_generator_recursion(2,2,$dcf,$hb);

application "fulton";
$n = 51;
$q = 37;
$cqs = new CyclicQuotient(N=>$n,Q=>$q);
$dcf = $cqs->DUAL_CONTINUED_FRACTION;
$hb = $cqs->WEIGHT_CONE->HILBERT_BASIS;
@result = ();
for(my $i=1; $i<=$n; $i++){
   @vec = ();
   for(my $j=1; $j<=$n; $j++){
      my $test = module_generator_recursion($i,$j,$dcf,$hb);
      push @vec, $test->rows();
   }
   push @result, new Vector(@vec);
}
print new Matrix(@result);
print $cqs->EXT1_MATRIX == new Matrix(@result);


