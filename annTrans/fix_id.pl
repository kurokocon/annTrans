#!/usr/bin/env perl


use strict;

open (my $file, "<", $ARGV[0]);

my $id = <$file>;
print $id;
if ($id =~ /ID=align_(\d+)/) {
	$id = $1;
}

$id++;
while (<$file>) {
	if ($_ =~ s/ID=align_\d+/ID=align_$id/) {
		print $_;

	}
	$id++;
}
