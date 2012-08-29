#!/usr/bin/perl -w

package CerevisiaeNetwork;
use Network;
@ISA = ("Network");

$| = 1; # Autoflush after each output.


# LOAD REQUIRED MODULES
use Storable qw(dclone);
use Gene;
use List::Util qw(min max);

# LOAD OPTIONAL MODULES
my $hasIntermine = 1;
eval "use Webservice::InterMine 0.9809 'yeastmine.yeastgenome.org/yeastmine'; 1" or $hasIntermine = 0;
my $hasGoPerl = 1;
eval "use GO::Parser; 1" or $hasGoPerl = 0;




#
#
#
# UPDATE FUNCTIONS
#
#
#
sub updateGenes {
	my ($self, $genus) = @_;
	
	$genus = "Saccharomyces" unless defined $genus;
	
	if(!$hasIntermine) {
		print "Cannot find Intermine 0.9809 Library. Resorting to file-based gene update.\n";
		print "Whoops... This feature is not implemented yet, sorry!\nPlease install Intermine.\n";
		return;
	} else {
		print "Querying Yeastmine database...\n";
	}
	
	# GENE PROPERTIES
	my $properties = Webservice::InterMine->new_query(class => 'Gene');
	$properties->add_view(qw/
					 primaryIdentifier
					 name
					 featureType
					 description
					 featAttribute
					 headline
					 length
					 ncbiGeneNumber
					 qualifier
					 secondaryIdentifier
					 sgdAlias
					 status
					 symbol
					 chromosome.primaryIdentifier
					 
					 locations.start
					 locations.end
					 locations.strand
					 
					 
					 organism.genus
					 /)->where('organism.genus' => $genus);
	# and filter the results for NULL values
	$properties->add_constraint("featureType","IS NOT NULL");
	$properties->add_constraint("primaryIdentifier","IS NOT NULL");
	$properties->add_constraint("chromosome.primaryIdentifier","!=","chrmt");
	$properties->add_constraint("featureType","!=","not in systematic sequence of S288C");
	
	
	# UPSTREAM INTERGENIC REGION
	my $upstream = Webservice::InterMine->new_query(class => 'Gene');
	# define upstream regions
	$upstream->add_view(qw/
					 primaryIdentifier
					 
					 upstreamIntergenicRegion.length
					 upstreamIntergenicRegion.chromosome.primaryIdentifier
					 upstreamIntergenicRegion.locations.start
					 upstreamIntergenicRegion.locations.end
					 
					 organism.genus
					 /)->where('organism.genus' => $genus);
	
	
	# DOWNSTREAM INTERGENIC REGION
	my $downstream = Webservice::InterMine->new_query(class => 'Gene');
	# define downstream regions
	$downstream->add_view(qw/
					 primaryIdentifier
					 
					 downstreamIntergenicRegion.locations.start
					 downstreamIntergenicRegion.locations.end
					 
					 organism.genus
					 /)->where('organism.genus' => $genus);
	
	
	
	# GENE ONTOLOGY TERMS
	my $go = Webservice::InterMine->new_query(class => 'Gene');
	# define downstream regions
	$go->add_view(qw/
					 primaryIdentifier
					 
					 goAnnotation.ontologyTerm.identifier
					 goAnnotation.ontologyTerm.name
					 
					 
					 organism.genus
					 /)->where('organism.genus' => $genus);
	
	
	
	# HOMOLOGUES
	my $homologues = Webservice::InterMine->new_query(class => 'Gene');
	# define downstream regions
	$homologues->add_view(qw/
					 primaryIdentifier
					 
					 homologues.homologue.primaryIdentifier
					 homologues.homologue.secondaryIdentifier
					 
					 
					 organism.genus
					 /)->where('organism.genus' => $genus);
	
	
	# INTERACTIONS
	my $interactions = Webservice::InterMine->new_query(class => 'Gene');
	# define downstream regions
	$interactions->add_view(qw/
					 primaryIdentifier
					 
					 interactions.interactingGenes.primaryIdentifier
					 interactions.interactionType
					 
					 
					 organism.genus
					 /)->where('organism.genus' => $genus);
	
	# PATHWAYS
	my $pathways = Webservice::InterMine->new_query(class => 'Gene');
	# define downstream regions
	$pathways->add_view(qw/
					 primaryIdentifier
					 
					 pathways.name
					 pathways.identifier
					 
					 organism.genus
					 /)->where('organism.genus' => $genus);
	
	# PUBLICATIONS
	my $publications = Webservice::InterMine->new_query(class => 'Gene');
	# define downstream regions
	$publications->add_view(qw/
					 primaryIdentifier
					 
					 publications.citation
					 publications.pubMedId
					 
					 organism.genus
					 /)->where('organism.genus' => $genus);
	
	my %translateRoman = (
		"I" => 1, "II" => 2, "III" => 3, "IV" => 4,
		"V" => 5, "VI" => 6, "VII" => 7, "VIII" => 8,
		"IX" => 9, "X" => 10, "XI" => 11, "XII" => 12,
		"XIII" => 13, "XIV" => 14, "XV" => 15, "XVI" => 16,
		"M" => "M"
	);
	
	
	
	print "\tMapping properties...\n";
	my $gene_ix = 0;
	my %genes;
	for (my $i = $properties->iterator(size => $properties->count); my $row = <$i>;) {
		my $gene = new Gene($gene_ix++,$row->{'primaryIdentifier'});
		
		$gene->set("secondaryIdentifier",$row->{'secondaryIdentifier'});
		$gene->set("alias",$row->{'sgdAlias'});
		$gene->set("ncbiGeneNumber",$row->{'ncbiGeneNumber'});
		$gene->set("name",$row->{'name'});
		$gene->set("type",$row->{'featureType'});
		$gene->set("description",$row->{'description'});
		$gene->set("attribute",$row->{'featAttribute'});
		$gene->set("headline",$row->{'headline'});
		$gene->set("length",$row->{'length'});
		$gene->set("qualifier",$row->{'qualifier'});
		$gene->set("status",$row->{'status'});
		$gene->set("symbol",$row->{'symbol'});
		my $chr = $row->{'chromosome.primaryIdentifier'};
		$chr =~ s/^chr//g;
		#$chr =~ s/^0//;
		$chr = $translateRoman{$chr};
		
		my $strand = $row->{'locations.strand'};
		$gene->set("chromosome",$chr);
		if($strand eq "-1") {
			$gene->set("tls",$row->{'locations.end'});
			$gene->set("tlt",$row->{'locations.start'});
		} else {
			$gene->set("tls",$row->{'locations.start'});
			$gene->set("tlt",$row->{'locations.end'});
		}
		$gene->set("start",$row->{'locations.start'});
		$gene->set("end",$row->{'locations.end'});
		$gene->set("strand",$strand);

		#$gene->set("sequence",$row->{'sequence.residues'});
		
		if($chr =~ /\d+/) {
			$genes{$gene->id} = $gene;
		} else {
			print $chr."\n";
			print "\t|--Excluding ".$gene->id." (not in nucleus)\n";
		}
		
	}
	
	print "\tMapping promoters...\n";
	for (my $i = $upstream->iterator(size => $upstream->count); my $row = <$i>;) {
		my $id = $row->{'primaryIdentifier'};
		
		my $gene = $genes{$id};
		if(defined $gene) {
			my $strand = $gene->get("strand");
			if($strand eq "-1") {
				$gene->set("promoter","start",$row->{'upstreamIntergenicRegion.locations.end'});
				$gene->set("promoter","end",$row->{'upstreamIntergenicRegion.locations.start'});
			} else {
				$gene->set("promoter","start",$row->{'upstreamIntergenicRegion.locations.start'});
				$gene->set("promoter","end",$row->{'upstreamIntergenicRegion.locations.end'});
			}
		}
	}
	
	print "\tMapping downstream regions...\n";
	for (my $i = $downstream->iterator(size => $downstream->count); my $row = <$i>;) {
		my $id = $row->{'primaryIdentifier'};
		
		my $gene = $genes{$id};
		if(defined $gene) {
			$gene->set("downstream","start",$row->{'downstreamIntergenicRegion.locations.start'});
			$gene->set("downstream","end",$row->{'downstreamIntergenicRegion.locations.end'});
		}
	}
	
	
	print "\tMapping Gene Ontology terms...\n";
	for (my $i = $go->iterator(size => $go->count); my $row = <$i>;) {
		my $id = $row->{'primaryIdentifier'};
		
		my $gene = $genes{$id};
		if(defined $gene) {
			#my $goId = $row->{'goAnnotation.ontologyTerm.identifier'};
			
			$gene->set("goTerms",$row->{'goAnnotation.ontologyTerm.identifier'},$row->{'goAnnotation.ontologyTerm.name'});
	
			
		}
	}
	
	print "\tMapping homologues...\n";
	for (my $i = $homologues->iterator(size => $homologues->count); my $row = <$i>;) {
		my $id = $row->{'primaryIdentifier'};
		
		my $gene = $genes{$id};
		if(defined $gene) {
			$gene->set("homologues",$row->{'homologues.homologue.primaryIdentifier'},$row->{'homologues.homologue.secondaryIdentifier'});
		}
	}
	
	print "\tMapping interactions...\n";
	for (my $i = $interactions->iterator(size => $interactions->count); my $row = <$i>;) {
		my $id = $row->{'primaryIdentifier'};
		
		my $gene = $genes{$id};
		if(defined $gene) {
			$gene->set('interactions',$row->{'interactions.interactingGenes.primaryIdentifier'},$row->{'interactions.interactionType'});
		}
	}
	
	print "\tMapping pathways...\n";
	for (my $i = $pathways->iterator(size => $pathways->count); my $row = <$i>;) {
		my $id = $row->{'primaryIdentifier'};
		
		my $gene = $genes{$id};
		if(defined $gene) {
			$gene->set('pathways',$row->{'pathways.identifier'},$row->{'pathways.name'});
		}
	}
	
	print "\tMapping publications...\n";
	for (my $i = $publications->iterator(size => $publications->count); my $row = <$i>;) {
		my $id = $row->{'primaryIdentifier'};
		
		my $gene = $genes{$id};
		if(defined $gene) {
			if(defined $row->{'publications.pubMedId'}) {
				$gene->set('publications',$row->{'publications.citation'},$row->{'publications.pubMedId'});
			} else {
				$gene->set('publications',$row->{'publications.citation'},"NA");
			}
			
		}
	}
	
	print "\tFilling missing values ...\n";
	foreach my $gene (values %genes) {
		if(!defined $gene->get('goTerms')) {
			my %gos;
			$gene->set('goTerms',\%gos);
		}
		if(!defined $gene->get('homologs')) {
			my %homologues;
			$gene->set('homologues',\%homologues);
		}
		if(!defined $gene->get('interactions')) {
			my %interactions;
			$gene->set('interactions',\%interactions);
		}
		if(!defined $gene->get('pathways')) {
			my %pathways;
			$gene->set('pathways',\%pathways);
		}
		if(!defined $gene->get('publications')) {
			my %publications;
			$gene->set('publications',\%publications);
		}
		#print $gene->to_string."\n"; 
	}
	
	#print scalar(keys %genes)."\n";
	
	print "\tDone.\n";
	
	$self->genes(\%genes);
	$self->build_gene_map;
}



sub updateNodes {
	my ($self) = @_;
	
	
	my $counterForOutputOnly = 0;
	
	$self->addGeneAssignmentStrategy("orf");
	$self->addGeneAssignmentStrategy("promoter");
	$self->addGeneAssignmentStrategy("tls");
	$self->addGeneAssignmentStrategy("largestOrf");
	$self->addGeneAssignmentStrategy("largestPromoter");
	
	
	my @unmappable = ();
	foreach my $gene (values %{$self->genes}) {
		my $geneChromosome = $gene->chromosome;
		my $tls = $gene->tls;
		my $strand = $gene->get("strand");
		my ($geneStart,$geneEnd) = ($gene->start,$gene->end);
		my ($promoterStart,$promoterEnd) = ($gene->upstreamStart,$gene->upstreamEnd);
		
		#fill undefined promoter regions
		if(!defined $promoterStart || !defined $promoterEnd) {
			my $promoterStartCandidate = 0;
			$promoterStartCandidate = $self->chromosomeLength($geneChromosome) if $strand == -1;
			
			foreach my $neighbor (values %{$self->genes}) {
				if($neighbor->chromosome == $geneChromosome && $gene->id ne $neighbor->id) {
					my $distStart = $tls - $neighbor->start;
					my $distEnd= $tls - $neighbor->end;
					
					if($strand eq "-1") {
						next if ($distStart > 0 || $distEnd > 0);
						my $possibleStart = $neighbor->start;
						if($distStart < $distEnd) {
							$possibleStart = $neighbor->end;
						}
						
						if(abs($tls-$possibleStart) < abs($tls-$promoterStartCandidate)) {
							$promoterStartCandidate = $possibleStart;
						}
					} else {
						next if ($distStart < 0 || $distEnd < 0);
						my $possibleStart = $neighbor->start;
						if($distStart > $distEnd) {
							$possibleStart = $neighbor->end;
						}
						
						if(abs($tls-$possibleStart) < abs($tls-$promoterStartCandidate)) {
							$promoterStartCandidate = $possibleStart;
						}
					}
					
				}
			}
			
			$gene->set("promoter","start",$promoterStartCandidate);
			$gene->set("promoter","end",$tls-1) if $strand == 1;
			$gene->set("promoter","end",$tls+1) if $strand == -1;
			
			($promoterStart,$promoterEnd) = ($gene->upstreamStart,$gene->upstreamEnd);
			
			#print $gene->id.": ".$gene->promoterStart."->".$gene->promoterEnd."->".$tls."->".$gene->tlt." <----------\n";
		} else {
			#print $gene->id.": ".$gene->promoterStart."->".$gene->promoterEnd."->".$tls."->".$gene->tlt."\n";
		}
		
		my $largestOrfNode;
		my $largestOrfSize = 0;;
		my $largestPromoterNode;
		my $largestPromoterSize = 0;
		
		foreach my $node (@{$self->nodes}) {
			my $nodeChromosome = $node->chromosome;
			my ($nodeStart,$nodeEnd) = ($node->start,$node->end);
			
			
			# by orf
			if($geneChromosome == $nodeChromosome) {
				
				my $length = 0;
				my $orf = 0;
				
				# overlapping on left
				if($geneStart <= $nodeStart && $geneEnd >= $nodeStart && $geneEnd <= $nodeEnd) {
					$length = $geneEnd - $nodeStart;
					$orf = 1;
				# contained completely	
				} elsif ($geneStart >= $nodeStart && $geneStart <= $nodeEnd && $geneEnd >= $nodeStart && $geneEnd <= $nodeEnd) {
					$length = $geneEnd - $geneStart;
					$orf = 1;
				# overlapping on right
				} elsif ($geneStart >= $nodeStart && $geneStart <= $nodeEnd && $geneEnd >= $nodeEnd) {
					$length = $nodeEnd - $geneStart;
					$orf = 1;
				# overlapping left and right
				} elsif($geneStart <= $nodeStart && $geneEnd >= $nodeEnd) {
					$length = $nodeEnd - $nodeStart;
					$orf = 1;
				}
				
				if($orf) {
					$node->addGene("orf",$gene);
					#$gene->parentNode($node);
					if($length > $largestOrfSize) {
						$largestOrfNode = $node;
						$largestOrfSize = $length;
					}
				}
			}
			
			# by promoter
			if(defined $promoterStart && defined $promoterEnd) {
				if($geneChromosome == $nodeChromosome) {
					my $length = 0;
					my $promoter = 0;
					
					# overlapping on left
					if($promoterStart <= $nodeStart && $promoterEnd >= $nodeStart && $promoterEnd <= $nodeEnd) {
						$length = $promoterEnd - $nodeStart;
						$promoter = 1;
					# contained completely	
					} elsif ($promoterStart >= $nodeStart && $promoterStart <= $nodeEnd && $promoterEnd >= $nodeStart && $promoterEnd <= $nodeEnd) {
						$length = $promoterEnd - $promoterStart;
						$promoter = 1;
					# overlapping on right
					} elsif ($promoterStart >= $nodeStart && $promoterStart <= $nodeEnd && $promoterEnd >= $nodeEnd) {
						$length = $nodeEnd - $promoterStart;
						$promoter = 1;
					# overlapping left and right
					} elsif($promoterStart <= $nodeStart && $promoterEnd >= $nodeEnd) {
						$length = $nodeEnd - $nodeStart;
						$promoter = 1;
					}
					
					if($promoter) {
						$node->addGene("promoter",$gene);
						#$gene->parentNode($node);
						if($length > $largestPromoterSize) {
							$largestPromoterNode = $node;
							$largestPromoterSize = $length;
						}
					}

				}
			}
			
			# by TLS
			if($geneChromosome == $nodeChromosome &&
			   (($tls <= $nodeEnd && $tls >= $nodeStart))) {
				$node->addGene("tls",$gene);
				#$gene->parentNode($node);
			}
			
		}
		
		# assign largest orf
		if(defined $largestOrfNode) {
			$largestOrfNode->addGene("largestOrf",$gene);
			$gene->parentNode($largestOrfNode);
		} else {
			print $gene->to_string."\n";
		}
		
		# assign largest promoter
		if(defined $largestPromoterNode) {
			$largestPromoterNode->addGene("largestPromoter",$gene);
		} else {
			push(@unmappable,$gene);
		}
		
		
		$counterForOutputOnly++;
		print "\r";
		print "Mapping genes to nodes... ".int($counterForOutputOnly/scalar(keys %{$self->genes}) * 100)."%";
	}
	print "\n";
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




1;
