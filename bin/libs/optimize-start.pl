#!/usr/bin/perl
use strict;
use warnings;

use 5.010;
use CGI::Carp qw(fatalsToBrowser);
use lib './libs/';
use OptimizerMain;

#Launch application
my $webapp = OptimizerMain -> new();
$webapp->run();