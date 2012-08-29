#!/usr/bin/perl -w

package Node;

sub new {
	my $class = shift;
	my $self = {
		_ix => shift,
		_chromosome => shift,
		_start => shift,
		_end => shift,
		_mappable => shift,
		_midpoint => undef,
		
		_genes => {},
		_data => {},
		
		
		_neighbors => {},
		_codingRegions => [],
		
		_originalIx => undef,
		
		_sequence => undef,
		
		_defaultStrategy => "largestOrf",

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
sub get {
	my $self = shift;
	
	my $data = $self->{_data};
	while(@_ > 0) {
		my $key = shift;
		if(!exists $$data{$key}) {
			return undef;
		}
		$data = $$data{$key};
	}
	
	return $data;
}

# fill Node object with data
# last value is datum, previous are hash keys
sub set {
	my $self = shift;
	
	my $data = $self->{_data};
		
	my $key = shift;
	while(@_ > 1) {
		$$data{$key} = {} unless exists $$data{$key};
		
		$data = $$data{$key};
		$key = shift;
	}
	$$data{$key} = shift;
	
	return;
}


sub ix {
	my ( $self, $ix ) = @_;
	$self->{_ix} = $ix if defined $ix;
	return $self->{_ix};
}

sub id {
	my ( $self, $ix ) = @_;
	$self->{_ix} = $ix if defined $ix;
	return $self->{_ix};
}

sub originalIx {
	my ( $self, $originalIx ) = @_;
	$self->{_originalIx} = $originalIx if defined $originalIx;
	return $self->{_originalIx};
}

sub chromosome {
	my ( $self, $chromosome ) = @_;
	$self->{_chromosome} = $chromosome if defined $chromosome;
	return $self->{_chromosome};
}

sub midpoint {
	my ( $self, $midpoint ) = @_;
	$self->{_midpoint} = $midpoint if defined $midpoint;
	return $self->{_midpoint};
}

sub startEnd {
	my($self) = @_;
	return ($self->start,$self->end) if $self->start<$self->end;
	return ($self->end,$self->start);
}


sub start {
	my ( $self, $start ) = @_;
	$self->{_start} = $start if defined $start;
	return $self->{_start};
}


sub end {
	my ( $self, $end ) = @_;
	$self->{_end} = $end if defined $end;
	return $self->{_end};
}

sub isMappable {
	my ( $self, $isMappable ) = @_;
	$self->{_mappable} = $isMappable if defined $isMappable;
	return $self->{_mappable};
}

sub mappable {
	my ( $self, $isMappable ) = @_;
	$self->{_mappable} = $isMappable if defined $isMappable;
	return $self->{_mappable};
}

sub sequence {
	my ( $self, $sequence ) = @_;
	$self->{_sequence} = $sequence if defined $sequence;
	return $self->{_sequence};
}


sub genes {
	my ( $self, $strategy ) = @_;
	if(defined $strategy) {
		if(exists $self->{_genes}{$strategy}) {
			return $self->{_genes}{$strategy};
		} else {
			return {};
		}
		
	}
	return $self->{_genes};
}

sub neighbors {
	my ( $self, $neighbors ) = @_;
	$self->{_neighbors} = $neighbors if defined $neighbors;
	return $self->{_neighbors};
}


sub data {
	my ( $self, $name ) = @_;
	if(defined $name) {
		if(exists $self->{_data}{$name}) {
			return $self->{_data}{$name};
		} else {
			return undef;
		}
		
	}
	return $self->{_data};
}












#
#
#
# DATA MAPPING FUNCTIONS
#
#
#
sub gene_labels_to_node_labels {
	my ($self,$name,$cutoff) = @_;
	$name = "last" unless defined $name;
	
	
	my %genes = %{$self->genes};
	foreach my $strategy (keys %genes) {
		
		my %nodeData;
		foreach my $gene (values %{$genes{$strategy}}) {
			my $data = $gene->get($name);
			
			if($data) {
				if(ref($data) eq "HASH") {
					foreach my $k (keys %{$data}) {
						if($cutoff && $$data{$k} > $cutoff) {
							$nodeData{$name}{$k} = 0 unless exists $nodeData{$name}{$k};
							$nodeData{$name}{$k}++
						} else {
							$nodeData{$name}{$k} = 0 unless exists $nodeData{$name}{$k};
							$nodeData{$name}{$k}++
						}
					}
				} else {
					$nodeData{$name}{$k} = 0 unless exists $nodeData{$name}{$k};
					$nodeData{$name}{$k}++
				}
			}
		}
		$self->{_data}{$strategy} = \%nodeData;
	}
}

























sub codingRegions {
	my ( $self, $codingRegions ) = @_;
	$self->{_codingRegions} = $codingRegions if defined $codingRegions;
	
	my @codingRegions = @{$self->{_codingRegions}};
	
	
	#sort coding regions
	my $swapped = 0;
	do {
		$swapped = 0;
		for(my $i=1; $i<@codingRegions; $i++) {
			if($codingRegions[$i-1][0] > $codingRegions[$i][0]) {
				my @saved = @{$codingRegions[$i-1]};
				$codingRegions[$i-1] = $codingRegions[$i];
				$codingRegions[$i] = \@saved;
				$swapped = 1;
			} 
		}
	} while ($swapped);
	
	my @mergedCodingRegions = ();
	my $j = 0;
	while($j < scalar(@codingRegions)) {
		my $start = $codingRegions[$j][0];
		my $end = $codingRegions[$j][1];
		
		
		my $changed = 1;
		while($changed && $j<scalar(@codingRegions)-1) {
			$changed = 0;
			if($codingRegions[$j+1][0] >= $start && $codingRegions[$j+1][0] <= $end) {
				#either it is overlapping
				if($codingRegions[$j+1][1] >= $start && $codingRegions[$j+1][1] > $end) {
					$end = $codingRegions[$j+1][1];
					$j++;
				}
				#or it is contained entirely
				else {
					$j++;
				}
				$changed = 1;
			}
		}
		push(@mergedCodingRegions,[$start,$end]);
		$j++;
	}
	
	return \@mergedCodingRegions;
}


sub intergenicRegions {
	my ( $self ) = @_;
	my @codingRegions = @{ $self->codingRegions };
	
	my($nodeStart,$nodeEnd) = $self->startEnd;
	
	
	my @intergenic = ();
	if(@codingRegions > 0) {
		my $start = $nodeStart;
		my $end;
		my $j = 0;
		if($codingRegions[0][0] == $nodeStart) {
			$start = $codingRegions[0][1];
			$j++;
		}
		while($j < scalar(@codingRegions)) {
			push(@intergenic,[$start+1,$codingRegions[$j][0]-1]);
			$start = $codingRegions[$j][1];
			$j++;
		}
		if($start != $nodeEnd) {
			push(@intergenic,[$start+1,$nodeEnd-1]);
		}		
	} else {
		my @empty = ();
		return \@empty;
	}
	
	return \@intergenic;
}

sub addGene {
	my ( $self, $strategy, $gene ) = @_;
	$self->{_genes}{$strategy}{$gene->id} = $gene;
}


sub addCodingRegion {
	my ( $self, $codingRegion ) = @_;
	push(@{ $self->{_codingRegions} }, $codingRegion) if defined $codingRegion;
}


sub addNeighbor {
	my ( $self, $neighbor ) = @_;
	$self->{_neighbors}{$neighbor->ix} = $neighbor if defined $neighbor;
}


sub length {
	my ( $self ) = @_;
	my $length = 0;
	($self->{_start} < $self->{_end}) ? ($length = $self->{_end}-$self->{_start}) : ($length = $self->{_start}-$self->{_end});
	return $length;
}

sub degree {
	my ( $self ) = @_;
	my $degree = scalar(keys %{$self->{_neighbors} });
	return $degree;
}


1;
