application "fulton";

$c = new Cone(INPUT_RAYS=>[[1,0,1],[1,1,0],[1,2,0],[1,1,1],[1,0,2]]);
$dc = new Cone(INPUT_RAYS=>$c->FACETS);
$tv = new NormalToricVariety($dc);
print $tv->ISOLATED; 

$c = new Cone(INPUT_RAYS=>[[1,0,0],[1,0,-1],[1,1,-4],[1,3,-8]]);
$dc = new Cone(INPUT_RAYS=>$c->FACETS);
$tv = new NormalToricVariety($dc);
print $tv->ISOLATED; 

$c = new Cone(INPUT_RAYS=>[]);
$dc = new Cone(INPUT_RAYS=>$c->FACETS);
$tv = new NormalToricVariety($dc);
print $tv->ISOLATED; 
