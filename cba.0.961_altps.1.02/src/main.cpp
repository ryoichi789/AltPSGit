/*
 * main.cpp: used by all cba applications
 */
#include "app.h"
extern Cba::App* app;

int main(int argc, char *argv[])
{
	app->set_param_default();
	app->set_param(argv);
	return app->run();
}
