use Node;
use Gene;
use Edge;
use Network;

use Storable;
use List::Util qw(shuffle);

$| = 1; # Autoflush after each output.
my $PID = $$;

my $error = <<END;

Usage:

Arg 1: number of chromosomes
Arg 2: number of nodes
Arg 3: beta (rewiring probability)
Arg 4: average degree
Arg 5: output

END

die "$error" unless scalar(@ARGV) >= 0;




my $nChromosomes = $ARGV[0];
my $nNodes = $ARGV[1];
my $beta = $ARGV[2];
my $avgDegree = $ARGV[3];
my $output = $ARGV[4];


# BASIC node definitions
my @nodes = ();
my %nodes;
my %availableNodes;
my $nodesPerChromosome = int($nNodes/$nChromosomes);
my $chromosome = 0;
for (my $ix=0; $ix<$nNodes; $ix++) {
	$chromosome++ if $ix % $nodesPerChromosome == 0 && $nNodes-$ix >= $nodesPerChromosome;
	
	my $n = new Node($ix,$chromosome,0,100,1);
	push(@nodes,$n);
	
	$nodes{$chromosome} = [] unless defined $nodes{$chromosome};
	push(@{$nodes{$chromosome}},$n);
	$availableNodes{$n->ix} = $n;
}



my @neighborEdges = ();
my %neighborEdges;

# NEIGHBOR edges
foreach my $chromosome (keys %nodes) {
	my @chrNodes = @{$nodes{$chromosome}};
		
	for(my $i=0; $i<@chrNodes-1; $i++) {
		my $n1 = $chrNodes[$i];
		my $n2 = $chrNodes[$i+1];
		
		$neighborEdges{$n1->ix}{$n2->ix} = scalar(@neighborEdges);
		push(@neighborEdges,[$n1->ix,$n2->ix]);
		$neighborEdges{$n2->ix}{$n1->ix} = scalar(@neighborEdges);
		push(@neighborEdges,[$n2->ix,$n1->ix]);
	}
}




# INTRA chromosomal edges (Watts and Strogatz model)
my @allIntraEdges = ();
my %allIntraEdges;
foreach my $chromosome (keys %nodes) {
	my @intraEdges = ();
	my %intraEdges;
	my @chrNodes = shuffle(@{$nodes{$chromosome}});
	
	for(my $i=0; $i<@chrNodes; $i++) {
		my $n1 = $chrNodes[$i];
		
		for(my $j=$i+1; ($j < @chrNodes) && ($j <= $i+$avgDegree/2); $j++) {
			my $n2 = $chrNodes[$j];
			
			if(!$neighborEdges{$n1->ix}{$n2->ix}) {
				$intraEdges{$n1->ix}{$n2->ix} = scalar(@intraEdges);
				push(@intraEdges,[$n1->ix,$n2->ix]);
				$intraEdges{$n2->ix}{$n1->ix} = scalar(@intraEdges);
				push(@intraEdges,[$n2->ix,$n1->ix]);
			}
		}
	}
	
	#rewiring
	for(my $i=0; $i<@intraEdges; $i+=2) {
		my $ix1 = $intraEdges{$intraEdges[$i][0]}{$intraEdges[$i][1]};
		my $ix2 = $intraEdges{$intraEdges[$i][1]}{$intraEdges[$i][0]};
				
		if(rand() < 0.5) {
			my $tmp = $ix1;
			$ix1 = $ix2;
			$ix2 = $tmp;
		}
		
		my $n1 = $intraEdges[$ix1][0];
		
		if(rand() < $beta) {
			my $n3;
			
			my $decline = 1;
			while($decline) {
				my $n = $chrNodes[int(rand(scalar(@chrNodes)))];
				$n3 = $n->ix;
								
				if(!$neighborEdges{$n1}{$n3} && !$intraEdges{$n1}{$n3} && $n1 != $n3) {
					$decline = 0;
				}
			}
			
			$intraEdges{$n1}{$n3} = $intraEdges{$intraEdges[$i][0]}{$intraEdges[$i][1]};
			$intraEdges{$n3}{$n1} = $intraEdges{$intraEdges[$i][1]}{$intraEdges[$i][0]};
			delete $intraEdges{$intraEdges[$i][0]}{$intraEdges[$i][1]};
			delete $intraEdges{$intraEdges[$i][1]}{$intraEdges[$i][0]};

			$intraEdges[$ix1][0] = $n1;
			$intraEdges[$ix1][1] = $n3;
			$intraEdges[$ix2][0] = $n3;
			$intraEdges[$ix2][1] = $n1;
		}
	}
	
	push(@allIntraEdges,@intraEdges);
}


#INTER edges (Watts and Strogatz model)
my @interNodeSequence = ();

my @k = shuffle(keys %availableNodes);

# get a feasible node order first
my @chrs = (-1) x int($avgDegree/2);
while(@k > 0) {
	my $n;
	
	my $i = 0;
	my $decline = 1;
	
	while($i < @k && $decline) {
		$n = $availableNodes{$k[$i]};
		
		$decline = 0;
		for(my $j=0; $j<@chrs; $j++) { $decline = 1 if $n->chromosome == $chrs[$j]; }
		$i++;
		
	}
	
	
	if($decline) {
		my $n = $availableNodes{$k[0]};
		splice(@k,0,1);
	} else {
		splice(@k,$i-1,1);
	}
	
	for(my $j=0; $j<@chrs-1; $j++) { $chrs[$j] = $chrs[$j+1]; }
	$chrs[scalar(@chrs)-1] = $n->chromosome;
	
	push(@interNodeSequence,$n);
	
}


# then start wiring and rewiring
my @interEdges = ();
my %interEdges;

for(my $i=0; $i<@interNodeSequence; $i++) {
	my $n1 = $interNodeSequence[$i];
	my $chr1 = $n1->chromosome;
	
	for(my $j=$i+1; ($j < @interNodeSequence) && ($j <= $i+$avgDegree/2); $j++) {
		my $n2 = $interNodeSequence[$j];
		my $chr2 = $n2->chromosome;
		
		if($chr1 != $chr2 && !$interEdges{$n1->ix}{$n2->ix}) {
			$interEdges{$n1->ix}{$n2->ix} = scalar(@interEdges);
			push(@interEdges,[$n1->ix,$n2->ix]);
			$interEdges{$n2->ix}{$n1->ix} = scalar(@interEdges);
			push(@interEdges,[$n2->ix,$n1->ix]);
		}
	}
}

#rewiring
for(my $i=0; $i<@interEdges; $i+=2) {
	my $ix1 = $interEdges{$interEdges[$i][0]}{$interEdges[$i][1]};
	my $ix2 = $interEdges{$interEdges[$i][1]}{$interEdges[$i][0]};
			
	if(rand() < 0.5) {
		my $tmp = $ix1;
		$ix1 = $ix2;
		$ix2 = $tmp;
	}
	
	my $n1 = $interEdges[$ix1][0];
	my $n1n = $availableNodes{$n1};
	
	if(rand() < $beta) {
		my $n3;
		
		my $decline = 1;
		while($decline) {
			my $n = $interNodeSequence[int(rand(scalar(@interNodeSequence)))];
			$n3 = $n->ix;
			
			
			if($n1 != $n3 && $n1n->chromosome != $n->chromosome) {
				$decline = 0;
			}
		}
		
		$interEdges{$n1}{$n3} = $interEdges{$interEdges[$i][0]}{$interEdges[$i][1]};
		$interEdges{$n3}{$n1} = $interEdges{$interEdges[$i][1]}{$interEdges[$i][0]};
		delete $interEdges{$interEdges[$i][0]}{$interEdges[$i][1]};
		delete $interEdges{$interEdges[$i][1]}{$interEdges[$i][0]};

		$interEdges[$ix1][0] = $n1;
		$interEdges[$ix1][1] = $n3;
		$interEdges[$ix2][0] = $n3;
		$interEdges[$ix2][1] = $n1;
	}
}





my @edges = ();
my $edgeIx = 0;


open(OUT,">$output.edgeList.txt") or die $!;

for(my $i=0; $i<@neighborEdges; $i+=2) {
	print OUT $neighborEdges[$i][0]."\t".$neighborEdges[$i][1]."\tneighbor\n";
	
	my $n1 = $availableNodes{$neighborEdges[$i][0]};
	my $n2 = $availableNodes{$neighborEdges[$i][1]};
	my $e = new Edge($edgeIx++,$n1,$n2,100,1,1,0);
	$e->type("neighbor");
	push(@edges,$e);
}

for(my $i=0; $i<@allIntraEdges; $i+=2) {
	print OUT $allIntraEdges[$i][0]."\t".$allIntraEdges[$i][1]."\tintra\n";
	
	my $n1 = $availableNodes{$allIntraEdges[$i][0]};
	my $n2 = $availableNodes{$allIntraEdges[$i][1]};
	my $e = new Edge($edgeIx++,$n1,$n2,100,1,1,0);
	$e->type("intra");
	push(@edges,$e);
}

for(my $i=0; $i<@interEdges; $i+=2) {
	print OUT $interEdges[$i][0]."\t".$interEdges[$i][1]."\tinter\n";
	
	my $n1 = $availableNodes{$interEdges[$i][0]};
	my $n2 = $availableNodes{$interEdges[$i][1]};
	my $e = new Edge($edgeIx++,$n1,$n2,100,1,1,0);
	$e->type("inter");
	push(@edges,$e);
}

close(OUT);




my $network = new Network(0.0001,\@nodes,\@edges);
store($network,"$output");

$network->corners;

exit;











