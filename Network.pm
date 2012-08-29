#!/usr/bin/perl -w

package Network;

$| = 1; # Autoflush after each output.


# LOAD REQUIRED MODULES
use Storable qw(dclone);
use Gene;
use List::Util qw(min max);




sub new {
	my $class = shift;
	my $self = {
		_fdr => shift,
		_nodes => shift,
		_edges => shift,
		
		_genes => {},
		
		_chromosomeLengths => [],
		_geneMap => {},
	};
	
	bless $self, $class;
}




#
#
#
# BASIC GETTERS / SETTERS
#
#
#
sub fdr {
	my ( $self, $fdr ) = @_;
	$self->{_fdr} = $fdr if defined $fdr;
	return $self->{_fdr};
}

sub nodes {
	my ( $self, $nodes ) = @_;
	$self->{_nodes} = $nodes if defined $nodes;
	return $self->{_nodes};
}

sub edges {
	my ( $self, $edges ) = @_;
	$self->{_edges} = $edges if defined $edges;
	return $self->{_edges};
}


sub genes {
	my ( $self, $genes ) = @_;
	$self->{_genes} = $genes if defined $genes;
	return $self->{_genes};
}

sub chromosomeLengths {
	my ( $self, $length ) = @_;
	$self->{_chromosomeLengths} = $length if defined $length;
	return $self->{_chromosomeLengths};
}

sub chromosomeLength {
	my ( $self, $chr ) = @_;
	my %lengths = %{$self->{_chromosomeLengths}};
	return $lengths{$chr};
}

sub edgeList {
	my ($self) = @_;
	
	my @edgeList = ();
	foreach my $e (@{$self->edges}) {
		my $n1 = $e->source;
		my $n2 = $e->sink;
		
		push(@edgeList, [$n1->ix,$n2->ix]);
	}
	
	return \@edgeList;
}





#
#
#
# ADVANCED GETTERS / SETTERS
#
#
#
sub nodeGroups {
	my ($self,$file,$oneGroupPerLine) = @_;
	$oneGroupPerLine = 1 unless defined $oneGroupPerLine;
	
	my @nodeGroups;
	my @names;
	
	
	if(defined $file) {
		open(IN,"$file") or die $!;
		
		my @group = ();
		my $groupNumber = 0;
		my $name = "0";
		my %names;
		while(<IN>) {
			chomp;
			if(!$oneGroupPerLine) {
				if($_ =~ /^#(.*)/) {
					push(@nodeGroups, \@group) if @group > 0;
					push(@names, $name) if @group > 0;
					@group = ();
					
					$name = $1;
					chomp($name);
					$name =~ s/\s//g;
					$name = $groupNumber++ if $name !~ /\w/;
					
					next;
				}
				
				$_ =~ s/,$//;
				$s =~ s/[^\d,]//g;
				@items = split(",",$s);
				
				push(@group,@items);
			} else {
				my $line = $_;
				if($line =~ /^#(.*)/) {
					$name = $1;
					chomp($name);
					$name =~ s/\s//g;
					next;
				}
				
				
				
				# make sure nodes are delimited by commas
				# enforce if not
				$line =~ s/\D/,/g;
				$line =~ s/(,)\1+/,/g;

				
				my @nodeIxs = split(",",$line);
				
				my @thisGroup = ();
				# add nodes from network
				foreach my $ix (@nodeIxs) {
					my $n = $self->getNode($ix);
					push(@thisGroup,$n) if defined $n;
				}
				push(@nodeGroups,\@thisGroup) if (@thisGroup > 0);
				
				if(exists $names{$name} || $name !~ /\w/) {
					$name = $groupNumber++
				}
				push(@names,$name);
				$names{$name} = 1;
			}
		}
	}
	return (\@nodeGroups, \@names);
}


sub getGene {
	my ($self,$gene) = @_;
	
	return $self->{_geneMap}{lc($gene)} if exists $self->{_geneMap}{lc($gene)};
	return undef;
}

sub getNode {
	my ($self,$node) = @_;
	
	my @nodes = @{$self->nodes};
	
	if($node =~ /^\d+$/) {
		return $nodes[$node];
	}
	return undef;
}


sub build_gene_map {
	my ($self) = @_;
	
	$self->{_geneMap} = {};
	
	foreach my $g (values %{$self->genes}) {
		$self->{_geneMap}{lc($g->id)} = $g;
		$self->{_geneMap}{lc($g->ix)} = $g;
		$self->{_geneMap}{lc($g->get("name"))} = $g if defined $g->get("name");
		$self->{_geneMap}{lc($g->get("secondaryIdentifier"))} = $g if defined $g->get("secondaryIdentifier");
		$self->{_geneMap}{lc($g->get("alias"))} = $g if defined $g->get("alias");
		$self->{_geneMap}{lc($g->get("symbol"))} = $g if defined $g->get("symbol");
		
		if(defined $g->get("secondaryIdentifier")) {
			my $sec = $g->get("secondaryIdentifier");
			if($sec =~ /(Y\w\w\d\d\d\w)-(\w)/) {
				my $alt = "".$1.$2;
				$self->{_geneMap}{lc($alt)} = $g;
			}
		}
	}
	
	return $self->{_geneMap};
}



sub adjacencyList {
	my ($self) = @_;
	
	my %edges;
	my $counter = 0;
	foreach my $e (@{$self->edges}) {
		my $source = ($e->source)->ix;
		my $sink = ($e->sink)->ix;
		
		if($sink == 4 || $source == 4) {
			$counter++;
		}
		
		$edges{$source}{$sink} = 1 if $source < $sink;
		$edges{$sink}{$source} = 1 if $source > $sink;
	}
	
	my @adjacency = ();
	for(my $i=0; $i<@{$self->nodes}; $i++) {
		$adjacency[$i] = [];
		if(exists $edges{$i}) {
			foreach my $j (sort(keys %{$edges{$i}})) {
				push(@{$adjacency[$i]},$j);
			}
		}
	}
	
	return \@adjacency;
}


sub edgesHash {
	my ($self) = @_;
	
	my %edges;
	my $counter = 0;
	foreach my $e (@{$self->edges}) {
		my $source = ($e->source)->ix;
		my $sink = ($e->sink)->ix;
		
		$edges{$source}{$sink} = 1;
		$edges{$sink}{$source} = 1;
	}
	
	return \%edges;
}





# GETTING SPECIFIC SUBNETWORKS
sub split_centromeric {
	# cutoff is one-sided, i.e. $cutoff to left and right of centromere
	my ($self,$cutoff,$centromerePointer) = @_;
	my @centromeres = ();
	if($centromerePointer) {
		@centromeres = @$centromerePointer;
	} else {
		@centromeres = (151524,238265,114443,449766,152046,148567,496979,105645,355687,436366,440188,150888,268090,628817,326643,556015);
	}
	
	my @nodes = @{$self->nodes};
	
	my %cNodes;
	my %nNodes;
	foreach my $n (@nodes) {
		my $chr = $n->chromosome;
		my $start = $n->start;
		my $end = $n->end;
		
		my $n_copy = new Node($n->ix,$n->chromosome,$n->start,$n->end,$n->isMappable,$n->midpoint);
		$n_copy->genes($n->genes);
		$n_copy->data($n->data);
		
		if(abs($start-$centromeres[$chr-1]) < $cutoff || abs($end-$centromeres[$chr-1]) < $cutoff) {
			$cNodes{$n_copy->ix} = $n_copy;
		} else {
			$nNodes{$n_copy->ix} = $n_copy;
		}
	}
	
	my @edges = @{$self->edges};
	
	my @cEdges = ();
	my @nEdges = ();
	my $nEdgeCounter = 0;
	my $cEdgeCounter = 0;
	foreach my $e (@edges) {
		my $source = $e->source;
		my $sink = $e->sink;
		
		if(exists $cNodes{$source->ix} && exists $cNodes{$sink->ix}) {
			my $edge = new Edge($cEdgeCounter++,$cNodes{$source->ix},$cNodes{$sink->ix},$e->frequency,$e->qValue);
			push(@cEdges,$edge);
			($cNodes{$source->ix})->addNeighbor($cNodes{$sink->ix});
			($cNodes{$sink->ix})->addNeighbor($cNodes{$source->ix});
		} elsif(exists $nNodes{$source->ix} && exists $nNodes{$sink->ix}) {
			my $edge = new Edge($nEdgeCounter++,$nNodes{$source->ix},$nNodes{$sink->ix},$e->frequency,$e->qValue);
			push(@nEdges,$edge);
			($nNodes{$source->ix})->addNeighbor($nNodes{$sink->ix});
			($nNodes{$sink->ix})->addNeighbor($nNodes{$source->ix});
		}
	}
	
	my @cNodes = values %cNodes;
	my @nNodes = values %nNodes;
	
	my $cNetwork = new Network($self->fdr,\@cNodes,\@cEdges);
	my $nNetwork = new Network($self->fdr,\@nNodes,\@nEdges);
	
	$cNetwork->genes($self->genes);
	$nNetwork->genes($self->genes);
	
	$cNetwork->geneAssignmentStrategies($self->geneAssignmentStrategies);
	$nNetwork->geneAssignmentStrategies($self->geneAssignmentStrategies);
	
	$cNetwork->chromosomeLengths($self->chromosomeLengths);
	$nNetwork->chromosomeLengths($self->chromosomeLengths);
	
	return($cNetwork,$nNetwork);
}


sub writeEdgeList {
	my ($self,$file) = @_;
	
	open(OUT,">$file") or die $!;
	
	my @edges = @{$self->edges};
	foreach my $e (@edges) {
		print OUT ($e->source)->ix."\t".($e->sink)->ix."\t".($e->weight)."\n";
	}
	
	close(OUT);
}










#
# function checked with igraph, works beautifully and fast
#
sub corners {
	my ($self,$verbose) = @_;
	$verbose = 0 unless defined $verbose;
	
	my @nodes = @{$self->nodes};
	my @edges = @{$self->edges};
	
	my $nTriangles = 0;
	
	#create edge list for instant lookup
	my %edges;
	my %types;
	foreach my $e (@edges) {
		my $source = $e->source;
		my $sink = $e->sink;
		$edges{$source->ix}{$sink->ix} = $e->type;
		$edges{$sink->ix}{$source->ix} = $e->type;
		$types{$e->type} = 1;
	}
	
	my @types = sort {$a cmp $b} keys %types;
	
	my %corners;
	for(my $i=0; $i<@types; $i++) {
		for(my $j=$i; $j<@types; $j++) {
			$corners{$types[$i]}{$types[$j]} = [];
			print $types[$i]." ".$types[$j]."\n" if $verbose;
		}
	}
	
	
	
	#search all triangles around node n
	for(my $i=0; $i<@nodes; $i++) {
		my $n1 = $nodes[$i];
		my $ix1 = $n1->ix;
		
		my @ixs2 = keys %{$edges{$ix1}};
		
		#speed up by only using direct neighbors
		foreach my $ix2 (@ixs2) {
			my @ixs3 = keys %{$edges{$ix2}};
			
			foreach my $ix3 (@ixs3) {
				
				#triangle found
				if($edges{$ix1}{$ix3}) {
					$nTriangles++;
					
					
					my $type1 = $edges{$ix1}{$ix2}; # 1. ix1--ix2
					my $type2 = $edges{$ix1}{$ix3}; # 2. ix1--ix3
					my $type3 = $edges{$ix2}{$ix3}; # 3. ix2--ix3
					
					# add 3 corners
					# pattern:
					# cornerstone, type1, type2
					
					# corner 1
					if($type1 le $type2) {
						push(@{$corners{$type1}{$type2}},[$ix1,$ix2,$ix3]);
					} else {
						push(@{$corners{$type2}{$type1}},[$ix1,$ix3,$ix2]);
					}
					
					# corner 2
					if($type1 le $type3) {
						push(@{$corners{$type1}{$type3}},[$ix2,$ix1,$ix3]);
					} else {
						push(@{$corners{$type3}{$type1}},[$ix2,$ix3,$ix1]);
					}
					
					# corner 3
					if($type2 le $type3) {
						push(@{$corners{$type2}{$type3}},[$ix3,$ix1,$ix2]);
					} else {
						push(@{$corners{$type3}{$type2}},[$ix3,$ix2,$ix1]);
					}
				}
			}
			
			#delete to prevent double counting of triangles
			delete $edges{$ix1}{$ix2};
			delete $edges{$ix2}{$ix1};
		}
	}
	
	if($verbose) {
		foreach my $type1 (keys %corners) {
			foreach my $type2 (keys %{$corners{$type1}}) {
				my @corners = @{$corners{$type1}{$type2}};
				print "------------ $type1 - $type2 ------------- ".scalar(@corners)."\n";
				foreach my $corner (@corners) {
					print $$corner[0]."-".$$corner[1]."-".$$corner[2]."\n";
				}
			}
		}
		print $nTriangles."\n";
	}
	
	return \%corners;
}



1;
