#!/usr/bin/perl
# 2006-05-11	Yo Matsuo

use strict;
$| = 1;

# parameters
my $FN = "";
my $NS = "Cba";
my $INC_GUARD = "";
my $CLASS_NAME = "";
my $NS_TMPL = "Cba";
my $INC_GUARD_TMPL = "CBA_CLASS_TMPL_H";
my $CLASS_NAME_TMPL = "ClassTmpl";
my $FTMPLH = "class_tmpl.h"; # "app_tmpl.h" for an application class
my $FTMPLCPP = "class_tmpl.cpp"; # "app_tmpl.cpp" for an application class

main();
1;

# main routine
sub main
{
	get_opt();
	set_param();
	make_src_from($FTMPLH);
	make_src_from($FTMPLCPP);
}

# usage (and exit)
sub usage()
{
	printf(STDERR "Usage: mksrc.pl [options]\n");
	printf(STDERR "\t-o: filename (lowercase letters + _; .h and .cpp files are created) [mandatory]\n");
	printf(STDERR "\t-s: namespace [Cba by default]\n");
	printf(STDERR "\t-a: y if making an application class [y/n; n by default]\n");
	exit(1);
}

# get command line options
sub get_opt
{
	my ($key, $value);
	while ($key = shift @ARGV) {
		$value = shift @ARGV || usage();
		if ($key eq "-o") {
			$FN = $value;
		} elsif ($key eq "-s") {
			$NS = $value;
		} elsif ($key eq "-a") {
			if (substr($value,0,1) eq "y") {
				$INC_GUARD_TMPL = "CBA_APP_TMPL_H";
				$CLASS_NAME_TMPL = "AppTmpl";
				$FTMPLH = "app_tmpl.h";
				$FTMPLCPP = "app_tmpl.cpp";
			}
		}
	}
	if (length($FN) == 0) {
		usage();
	}
}

# set parameters
sub set_param
{
	$INC_GUARD = uc($NS) . "_" . uc($FN) . "_H";
	my @words = split(/_/, $FN);
	for (@words) {
		$CLASS_NAME .= uc(substr($_,0,1)) . lc(substr($_,1));
	}
}

#
sub make_src_from
{
	my ($tmpl) = @_;
	my $ftmpl = "";
	my $fext = "";
	if ($tmpl =~ /^(.+)\.([^\.]+)$/) {
		$ftmpl = $1;
		$fext = $2;
	}
	return if (length($ftmpl) == 0);
	return if (length($fext) == 0);
	my $fn = $FN . "." . $fext;

	open(FI, "$tmpl");
	open(FO, ">$fn");
	while (<FI>) {
		next if (/^\s*\/\/\//);
		s/$ftmpl/$FN/g;
		s/$INC_GUARD_TMPL/$INC_GUARD/g; # replace include guards
		s/namespace Cba/namespace $NS/g; # replace namespace
		s/Cba::/$NS::/g;
		if ($NS ne "Cba") {
			s/ public App/ public Cba::App/g;
		}
		s/$CLASS_NAME_TMPL/$CLASS_NAME/g; # replace class names
		s/$ftmpl/$FN/g; # replace application name (filename)
		printf(FO "%s", $_);
	}
	close(FI);
	close(FO);
}
