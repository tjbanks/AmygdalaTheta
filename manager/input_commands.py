import click
import logging
import os

from clint.textui import puts, colored, indent
from .input_functions import (print_test)

@click.group('input')
@click.pass_context
def cli(ctx):

    pass


@click.group('generate', help='Generate input for an experiment')
@click.pass_context
def generate(ctx):

    pass

@generate.command('all',help="Generate input for all experiments")
@click.pass_context
def input_generate_all(ctx):
    print_test()

@generate.command('exp0',help="Generate input for experiment 1")
@click.pass_context
def input_generate_exp0(ctx):
    print_test()

cli.add_command(generate)