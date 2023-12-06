package My_xgrm;
use Mojo::Base 'Mojolicious';

# This method will run once at server start
sub startup {
	my $self = shift;
	
	$ENV{MOJO_REVERSE_PROXY} = 1;
	$self->config(
		hypnotoad => {
			listen  => ['http://*:3010'],
			workers => 8,
			keep_alive_timeout => 300,
			websocket_timeout => 600,
			proxy => 1
		}
	);
	
	# Documentation browser under "/perldoc"
  	$self->plugin('PODRenderer');
  	#$self->plugin('PODViewer');

  	# Router
  	my $r = $self->routes;
	
	# Template names are expected to follow the template.format.handler scheme, with template defaulting to controller/action or the route name, format defaulting to html and handler to ep
	
  	# Normal route to controller
  	## Home
  	$r->get('/')->to(template=>'index', controller=>'action', action=>'index');
  	
  	## XGRm (and typos)
	$r->get('/XGRm')->to(template=>'index', controller=>'action', action=>'index');
	$r->get('/XGRM')->to(template=>'index', controller=>'action', action=>'index');
	$r->get('/xgrm')->to(template=>'index', controller=>'action', action=>'index');
	$r->get('/XGRmod')->to(template=>'index', controller=>'action', action=>'index');
	$r->get('/xgrmod')->to(template=>'index', controller=>'action', action=>'index');
	$r->get('/XGRMOD')->to(template=>'index', controller=>'action', action=>'index');
	$r->get('/XGRmodel')->to(template=>'index', controller=>'action', action=>'index');
	$r->get('/xgrmodel')->to(template=>'index', controller=>'action', action=>'index');
	$r->get('/XGRMODEL')->to(template=>'index', controller=>'action', action=>'index');

  	## howto
  	$r->get('/XGRm/howto')->to(template=>'XGRm_howto', format=>'html', handler=>'ep', controller=>'action', action=>'index');

  	#############################################
  	### GElyser
	$r->get('/XGRm/GElyser')->to(template=>'GElyser', format=>'html', handler=>'ep', controller=>'action', action=>'GElyser', post_flag=>0);
	$r->post('/XGRm/GElyser')->to(template=>'GElyser', format=>'html', handler=>'ep', controller=>'action', action=>'GElyser', post_flag=>1);
	
  	### GSlyser
	$r->get('/XGRm/GSlyser')->to(template=>'GSlyser', format=>'html', handler=>'ep', controller=>'action', action=>'GSlyser', post_flag=>0);
	$r->post('/XGRm/GSlyser')->to(template=>'GSlyser', format=>'html', handler=>'ep', controller=>'action', action=>'GSlyser', post_flag=>1);
	
  	### RElyser
	$r->get('/XGRm/RElyser')->to(template=>'RElyser', format=>'html', handler=>'ep', controller=>'action', action=>'RElyser', post_flag=>0);
	$r->post('/XGRm/RElyser')->to(template=>'RElyser', format=>'html', handler=>'ep', controller=>'action', action=>'RElyser', post_flag=>1);
	
  	### RSlyser
	$r->get('/XGRm/RSlyser')->to(template=>'RSlyser', format=>'html', handler=>'ep', controller=>'action', action=>'RSlyser', post_flag=>0);
	$r->post('/XGRm/RSlyser')->to(template=>'RSlyser', format=>'html', handler=>'ep', controller=>'action', action=>'RSlyser', post_flag=>1);
	
}

1;

