my $c = cube(3);
compare_data("cube3", dcch($c->VERTICES, unit_vector($c->CONE_DIM, 0), 2));

my $s = cube(2);
compare_data("cube2", dcch($s->VERTICES, unit_vector($s->CONE_DIM, 0), 2));

my $s12 = simplex(12);
compare_data("simplex12", dcch($s12->VERTICES, unit_vector($s12->CONE_DIM, 0), 11));
