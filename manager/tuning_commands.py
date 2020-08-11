import click
import logging
import os

from clint.textui import puts, colored, indent
from .input_functions import (print_test)

@click.group('tune')
@click.pass_context
def cli(ctx):
    pass

@click.group('clamp')
@click.pass_context
def cli_clamp(ctx):
    pass

@click.group('fi')
@click.pass_context
def cli_fi(ctx):
    pass

@click.group('syn')
@click.pass_context
def cli_syn(ctx):
    pass

cli.add_command(cli_clamp)
cli.add_command(cli_fi)
cli.add_command(cli_syn)

#######################################################################
## Individual Commands                                               ##
#######################################################################

def execute(command):
    print("Executing '" + command + "'")
    os.system(command)


#######################################################################
## Current Clamps                                                    ##
#######################################################################

@cli_clamp.command('Cell_A',help="Current Clamp Cell A")
@click.pass_context
def cli_clamp_cell_a(ctx):
    execute('nrngui ./tuning/current_clamp/Cell_A.hoc')




#######################################################################
## FI Curves                                                         ##
#######################################################################

@cli_fi.command('Cell_A',help="FI Curve for Cell A")
@click.pass_context
def cli_fi_cell_a(ctx):
    execute('bmtool util cell --template Cell_A fi')
  
@cli_fi.command('Cell_C',help="FI Curve for Cell C")
@click.pass_context
def cli_fi_cell_c(ctx):
    execute('bmtool util cell --template Cell_C fi')

@cli_fi.command('PV',help="FI Curve for PV Cell")
@click.pass_context
def cli_fi_pv(ctx):
    execute('bmtool util cell --template PV fi')

@cli_fi.command('CR',help="FI Curve for CR Cell")
@click.pass_context
def cli_fi_cr(ctx):
    execute('bmtool util cell --template CR fi')

@cli_fi.command('SOM',help="FI Curve for SOM Cell")
@click.pass_context
def cli_fi_som(ctx):
    execute('bmtool util cell --template SOM fi')


#######################################################################
## Synapse Testing                                                   ##
#######################################################################

@cli_syn.command('exp2syn_Cell_A',help="Exp2Syn Cell A")
@click.pass_context
def cli_exp2syn_cell_a(ctx):
    execute('nrngui ./tuning/synapses/exp2syn_Cell_A.hoc')
    

