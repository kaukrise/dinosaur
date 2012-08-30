#!/usr/bin/perl -w

package RandomNetwork;
use Network;
use Edge;
use Node;
use List::Util qw(sum shuffle);
use Time::HiRes qw(gettimeofday tv_interval);

$| = 1; # Autoflush after each output.





# algorithm is 2-step filtering based:
# 1. corner1 ---filter--> corner2
# 2. corner2 ---filter--> corner1
# thereby failures are kept to a minimun
# (or avoided completely if done well)
#
#
# swapping restrictions
#
# intra // intra:
#     chromosomes of all nodes have to be identical
# intra // inter:
#     chromosomes of root nodes have to be identical
# intra // neighbor
#     ---
# intra // intra-inter
#     chromosome of intra root node has to be identical
#     to either of the intra-inter leaf nodes
# intra // intra-neighbor or intra // inter-neighbor
#     ---
#
# inter // inter
#     chromosome of first corner root node has to be
#     different from both chromosomes of second corner
#     leaf nodes and vice versa
sub randomize {
	my ($network) = @_;
	
	my @edges = @{$network->edges};
	my %edges;
	foreach my $e (@edges) {
		my $source = $e->source;
		my $sink = $e->sink;
		$edges{$source->ix}{$sink->ix} = $e->type;
		$edges{$sink->ix}{$source->ix} = $e->type;
	}
	
	my @nodes = @{$network->nodes};
	my %chromosomes;
	foreach my $n (@nodes) {
		$chromosomes{$n->ix} = $n->chromosome;
	}
	
	my @corners = @{$network->corners};
	my %corners = %{ &assign_corner_types($network,\@corners)};
	
	# restructure corners by chromosome and
	# identify edges that don't participate in
	# triangles of any kind
	#my %corners_by_type = %{$network->corners};
	#my %corners;
	#my %counts;
	#foreach my $type1 (keys %corners_by_type) {
	#	foreach my $type2 (keys %{$corners_by_type{$type1}}) {
	#
	#		my @corners = @{$corners_by_type{$type1}{$type2}};
	#		$counts{$type1}{$type2} = scalar(@corners);
	#		foreach my $corner (@corners) {
	#			my $root = $nodes{$$corner[0]};
	#			my $leaf1 = $nodes{$$corner[1]};
	#			my $leaf2 = $nodes{$$corner[2]};
	#			
	#			$edges{$root->ix}{$leaf1->ix} = 0;
	#			$edges{$leaf1->ix}{$root->ix} = 0;
	#			$edges{$root->ix}{$leaf2->ix} = 0;
	#			$edges{$leaf2->ix}{$root->ix} = 0;
	#			
	#			$corners{$type1}{$type2}{$root->chromosome}{$leaf1->chromosome}{$leaf2->chromosome} = [] unless exists $corners{$type1}{$type2}{$root->chromosome}{$leaf1->chromosome}{$leaf2->chromosome};
	#			push(@{$corners{$type1}{$type2}{$root->chromosome}{$leaf1->chromosome}{$leaf2->chromosome}},$corner);
	#		}
	#	}
	#}
	#
	#my %singles;
	#foreach my $e (@edges) {
	#	my $source = $e->source;
	#	my $sink = $e->sink;
	#	
	#	$singles{$source->chromosome}{$sink->chromosome} = [] unless exists $singles{$source->chromosome}{$sink->chromosome};
	#	push(@{$singles{$source->chromosome}{$sink->chromosome}},[$source->ix,$sink->ix]);
	#}
	
	# possible combinations
	# neighbor neighbor
	# inter neighbor
	# inter inter
	# inter intra
	# intra neighbor
	# intra intra
	 
	my $nCorners = scalar(@corners);
	my $swaps = $nCorners;
	while($swaps > 0) {
		my $corner1_ix = int(rand($nCorners));
		my @corner1 = @{$corners[$corner1_ix]};
		
		# determine corner type
		my $t1 = $edges{$corner1[0]}{$corner1[1]};
		my $t2 = $edges{$corner1[0]}{$corner1[2]};
		my $rootChr = $chromosomes{$corner1[0]};
		my $leaf1Chr = $chromosomes{$corner1[1]};
		my $leaf2Chr = $chromosomes{$corner1[2]};
		
		# reduce search space to viable solutions
		if($t1 eq "intra" && $t2 eq "intra") {
			my @possible_swap_corners = ();
			# all intra-intra with the same chromosomes
			push(@possible_swap_corners, keys %{$corners{"intra"}{"intra"}{$rootChr}{$rootChr}{$rootChr}});
			
			# all inter-inter with the same root chromosome
			foreach my $chr2 (keys %{$corners{"inter"}{"inter"}{$rootChr}}) {
				foreach my $chr3 (keys %{$corners{"inter"}{"inter"}{$rootChr}{$chr2}}) {
					push(@possible_swap_corners, keys %{$corners{"inter"}{"inter"}{$rootChr}{$chr2}{$chr3}});
				}
			}
			
			# all inter-intra with the same root chromosome
			foreach my $chr2 (keys %{$corners{"inter"}{"intra"}{$rootChr}}) {
				foreach my $chr3 (keys %{$corners{"inter"}{"intra"}{$rootChr}{$chr2}}) {
					push(@possible_swap_corners, keys %{$corners{"inter"}{"intra"}{$rootChr}{$chr2}{$chr3}});
				}
			}
						
			my $found_swap_partner = 0;
			while(!$found_swap_partner && scalar(@possible_swap_corners) > 0) {
				my $corner2_swap_ix = int(rand(scalar(@possible_swap_corners)));
				my $corner2_ix = $possible_swap_corners[$corner2_swap_ix];
				my @corner2 = @{$corners[$corner2_ix]};
								
				# determine type and position
				my $t1_2 = $edges{$corner2[0]}{$corner2[1]};
				my $t2_2 = $edges{$corner2[0]}{$corner2[2]};
				my $rootChr_2 = $chromosomes{$corner2[0]};
				my $leaf1Chr_2 = $chromosomes{$corner2[1]};
				my $leaf2Chr_2 = $chromosomes{$corner2[2]};
				
				if($corner2[0] != $corner1[1] && $corner2[0] != $corner1[2] &&
				   $corner1[0] != $corner2[1] && $corner1[0] != $corner2[2]) {
					$found_swap_partner = 1;
					
					delete $edges{$corner1[0]}{$corner1[1]};
					delete $edges{$corner1[1]}{$corner1[0]};
					delete $edges{$corner2[0]}{$corner2[1]};
					delete $edges{$corner2[1]}{$corner2[0]};
					
					#here is the swap
					my $tmp_swap = $corner1[0];
					$corners[$corner1_ix][0] = $corner2[0];
					$corners[$corner2_ix][0] = $tmp_swap;
					if(rand() < 0.5) {
						my $tmp = $corners[$corner1_ix][1];
						$corners[$corner1_ix][1] = $corners[$corner1_ix][2];
						$corners[$corner1_ix][2] = $corners[$corner1_ix][1];
					}
					
					$edges{$corner1[0]}{$corner1[1]} = $t1;
					$edges{$corner1[0]}{$corner1[2]} = $t2;
					$edges{$corner2[0]}{$corner2[1]} = $t1_2;
					$edges{$corner2[0]}{$corner2[2]} = $t2_2;
				} else {
					splice(@possible_swap_corners,$corner2_swap_ix,1);
				}
			}
		}
		
		$swaps--;
	}
}





sub assign_corner_types {
	my ($network,$corners) = @_;
	
	my @nodes = @{$network->nodes};
	my %nodes;
	foreach my $n (@nodes) {
		$nodes{$n->ix} = $n;
	}
	
	my @edges = @{$network->edges};
	my %edges;
	foreach my $e (@edges) {
		my $source = $e->source;
		my $sink = $e->sink;
		$edges{$source->ix}{$sink->ix} = $e->type;
		$edges{$sink->ix}{$source->ix} = $e->type;
	}	
	
	my %corners;
	for(my $i=0; $i<@$corners; $i++) {
		my $corner = $$corners[$i];
		my $root = $$corner[0];
		my $leaf1 = $$corner[1];
		my $leaf2 = $$corner[2];
		
		my $rootChr = ($nodes{$root})->chromosome;
		my $leaf1Chr = ($nodes{$leaf1})->chromosome;
		my $leaf2Chr = ($nodes{$leaf2})->chromosome;
		
		my $t1 = $edges{$root}{$leaf1};
		my $t2 = $edges{$root}{$leaf2};
		
		if($t1 le $t2) {
			#$corners{$t1}{$t2}{$rootChr}{$leaf1Chr}{$leaf2Chr} = [] unless exists $corners{$t1}{$t2}{$rootChr}{$leaf1Chr}{$leaf2Chr};
			#push(@{$corners{$t1}{$t2}{$rootChr}{$leaf1Chr}{$leaf2Chr}},[$root,$leaf1,$leaf2]);
			$corners{$t1}{$t2}{$rootChr}{$leaf1Chr}{$leaf2Chr}{$i} = 1;
		} else {
			#$corners{$t2}{$t1}{$rootChr}{$leaf2Chr}{$leaf1Chr} = [] unless exists $corners{$t2}{$t1}{$rootChr}{$leaf2Chr}{$leaf1Chr};
			#push(@{$corners{$t2}{$t1}{$rootChr}{$leaf2Chr}{$leaf1Chr}},[$root,$leaf2,$leaf1]);
			$corners{$t2}{$t1}{$rootChr}{$leaf2Chr}{$leaf1Chr}{$i} = 1;
		}
	}
	
	
	return \%corners;
}


1;
