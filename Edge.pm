#!/usr/bin/perl -w

package Edge;

sub new {
	my $class = shift;
	my $self = {
		_ix => shift,
		_source => shift,
		_sink => shift,
		_frequency => shift,
		_weight => shift,
		_probability => shift,
		_qValue => shift,
		
		_type => undef,
		
		_originalIx => undef,
	};
	
	bless $self, $class;
}



sub ix {
	my ( $self, $ix ) = @_;
	$self->{_ix} = $ix if defined $ix;
	return $self->{_ix};
}

sub originalIx {
	my ( $self, $originalIx ) = @_;
	$self->{_originalIx} = $originalIx if defined $originalIx;
	return $self->{_originalIx};
}

sub source {
	my ( $self, $source ) = @_;
	$self->{_source} = $source if defined $source;
	return $self->{_source};
}

sub sink {
	my ( $self, $sink ) = @_;
	$self->{_sink} = $sink if defined $sink;
	return $self->{_sink};
}

sub frequency {
	my ( $self, $frequency ) = @_;
	$self->{_frequency} = $frequency if defined $frequency;
	return $self->{_frequency};
}

sub probability {
	my ( $self, $probability ) = @_;
	$self->{_probability} = $probability if defined $probability;
	return $self->{_probability};
}

sub qValue {
	my ( $self, $qValue ) = @_;
	$self->{_qValue} = $qValue if defined $qValue;
	return $self->{_qValue};
}

sub weight {
	my ( $self, $weight ) = @_;
	$self->{_weight} = $weight if defined $weight;
	return $self->{_weight};
}

sub type {
	my ( $self, $type ) = @_;
	$self->{_type} = $type if defined $type;
	return $self->{_type};
}


sub clone {
	my ( $self ) = @_;
	my $e = new Edge($self->ix,$self->source,$self->sink,$self->frequency,$self->qValue);
	
	$e->originalIx($self->ix);
	return $e;
}

sub emptyClone {
	my ( $self ) = @_;
	my $e = new Edge($self->ix,undef,undef,$self->frequency,$self->qValue);
	
	$e->originalIx($self->ix);
	return $e;
}

1;
