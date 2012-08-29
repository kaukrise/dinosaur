#!/usr/bin/perl -w

package Gene;

sub new {
	my $class = shift;
	my $self = {
		_ix => shift,
		_id => shift,
		_data => {},
		
		_parentNode => undef,
		_originalIx => undef,
		
		_neighbors => {}
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

# fill Gene object with data
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

sub to_string {
    my ($self, $showValues, $string, $level, $hashref) = @_;
	$string = "Index: ".$self->ix."\nID: ".$self->id."\n" unless defined $string;
	$showValues = 1 unless defined $showValues;
	$level = 0 unless defined $level;
	$hashref = $self->{_data} unless defined $hashref;
	
	my @keys = sort(keys(%$hashref));
    foreach my $k (@keys) {
		my $v = $$hashref{$k};
		
		my $prefix = "";
		for(my $l=1; $l<$level; $l++) {
			$prefix .= "    "
		}
		$prefix .= "|-- " unless $level == 0;
		
        if (ref($v) eq 'HASH') {
			$string .= $prefix."".$k."\n";
            $string = to_string($v,$showValues,$string,$level+1,$v);
        } else {
			if(defined $v) {
				$string .= "$prefix".""."$k -> $v\n" if $showValues;
			} else {
				$string .= "$prefix".""."$k -> NA\n" if $showValues;
			}
            
        }
    }
	
	return($string);
}






sub ix {
	my ( $self, $ix ) = @_;
	$self->{_ix} = $ix if defined $ix;
	return $self->{_ix};
}

sub id {
	my ( $self, $id ) = @_;
	$self->{_id} = $id if defined $id;
	return $self->{_id};
}

sub originalIx {
	my ( $self, $originalIx ) = @_;
	$self->{_originalIx} = $originalIx if defined $originalIx;
	return $self->{_originalIx};
}


sub parentNode {
	my($self,$parentNode) = @_;
	$self->{_parentNode} = $parentNode if defined $parentNode;
	return $self->{_parentNode};
}

sub chromosome {
	my($self,$chr) = @_;
	$self->set("chromosome",$chr) if defined $chr;
	return $self->get("chromosome");
}

sub tls {
	my($self,$tls) = @_;
	$self->set("tls",$tls) if defined $tls;
	return $self->get("tls");
}

sub tlt {
	my($self,$tlt) = @_;
	$self->set("tlt",$tlt) if defined $tlt;
	return $self->get("tlt");
}

sub start {
	my($self) = @_;
	my $start = ($self->tls < $self->tlt) ? $self->tls : $self->tlt;
	return $start;
}

sub end {
	my($self) = @_;
	my $end = ($self->tls > $self->tlt) ? $self->tls : $self->tlt;
	return $end;
}

sub length {
	my ($self) = @_;
	return $self->end - $self->start;
}

sub promoterStart {
	my($self) = @_;
	if(!defined $self->get("promoter","start")) {
		return;
	}
	
	return $self->get("promoter","start");
}

sub promoterEnd {
	my($self) = @_;
	if(!defined $self->get("promoter","end")) {
		return;
	}
	
	return $self->get("promoter","end");
}

sub upstreamStart {
	my($self) = @_;
	if(!defined $self->get("promoter","start") || !defined $self->get("promoter","end")) {
		return;
	}
	if($self->get("promoter","start") < $self->get("promoter","end")) {
		return $self->get("promoter","start");
	}
	return $self->get("promoter","end");
}

sub upstreamEnd {
	my($self) = @_;
	if(!defined $self->get("promoter","start") || !defined $self->get("promoter","end")) {
		return;
	}
	if($self->get("promoter","start") < $self->get("promoter","end")) {
		return $self->get("promoter","end");
	}
	return $self->get("promoter","start");
}





sub neighbors {
	my ( $self, $neighbors ) = @_;
	$self->{_neighbors} = $neighbors if defined $neighbors;
	return $self->{_neighbors};
}


sub addNeighbor {
	my ( $self, $neighbor ) = @_;
	$self->{_neighbors}{$neighbor->id} = $neighbor if defined $neighbor;
}

sub degree {
	my ($self) = @_;
	return scalar(keys %{$self->neighbors});
}

1;
