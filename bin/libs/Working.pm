package Working;

# The Working package provides a simple
# mechanism to keep browsers from timing out on
# long-running processes when they would
# otherwise receive no output.
# 
# At the beginning of the long-running quiet
# code, call Working::start(n), where n is the
# number of seconds until/between messages to the
# browser. Then, when the process is completed
# and you are ready to start sending "real"
# output to the browser again, call
# Working::stop().

# Thanks to Todd Lewis for providing this code. 

# Sample code to demo this package is:
# 
#   Working::start( 3 );
#   $loop = 10;
#   while ( $loop-- )
#     {
#       print "($loop) sleeping\n";
#       sleep 1;
#       Working::stop() if $loop < 5;
#     }

my $working_time;

sub stop
  {
    alarm 0;
  }

sub working
  {
    $| = 1;  # Disable buffering
    print "<!-- working... -->\n";
    alarm $working_time;
  }

sub start($)
  {
    $working_time = shift;
    $SIG{ALRM} = \&working;
    alarm $working_time;
  }
1;