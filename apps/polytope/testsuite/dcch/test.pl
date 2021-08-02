my $c = cube(3);
compare_data("cube3", dcch($c->VERTICES, unit_vector($c->CONE_DIM, 0), 2));

my $s = cube(2);
compare_data("cube2", dcch($s->VERTICES, unit_vector($s->CONE_DIM, 0), 2));

my $s12 = simplex(12);
compare_data("simplex12", dcch($s12->VERTICES, unit_vector($s12->CONE_DIM, 0), 11));

my $pos_orth = new Cone(INPUT_RAYS=>[[1,0,0],[0,1,0],[0,0,1]]);
compare_data("pos_orth0", dcch($pos_orth->VERTICES, new Vector([1,1,1]), 2));
compare_data("pos_orth1", dcch($pos_orth->VERTICES, new Vector([10,1,1]), 2));
compare_data("pos_orth2", dcch($pos_orth->VERTICES, new Vector([3,17,8]), 2));
compare_data("pos_orth3", dcch($pos_orth->VERTICES, new Vector([300,2,50]), 2));
