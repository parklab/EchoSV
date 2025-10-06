#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@Time: 2025/10/01
@Author: Yuwei Zhang
@Contact: yuwei_zhang@hms.harvard.edu
'''

import echosv, argparse, sys
from echosv.generate_chain import chain_main
from echosv.genotype import genotype_main
from echosv.match import match_main
from rich.console import Console

USAGE = f"""\
[bold]EchoSV v{echosv.__version__}[/] - Structural Variant Comparison and Merging Toolkit across Different Reference Genomes

[bold cyan]chain[/]  Generate liftover chain file between two reference genomes 
[bold cyan]genotype[/]  Collecting supporting reads of structural variants using long/short read BAMs
[bold cyan]match[/]  Compare two/multi VCF files and generate a comparison result
"""

class ArgumentParser(argparse.ArgumentParser):
    """
    Custom argument parser error
    """

    def error(self, message):
        """
        Check for similar commands before exiting
        """
        console = Console(stderr=True)
        console.print(f'{self.prog}: error: {message}')
        if message.startswith("argument CMD: invalid choice"):
            guess = echosv.help_unknown_cmd(sys.argv[1], ["chain", "genotype", "match"])
            if guess:
                console.print(f"\nThe most similar command is\n\t[bold][cyan]{guess}[/][/]\n")
        self.exit(2)

    def _print_message(self, message, file=None):
        """ pretty print """
        console = Console(stderr=True)
        console.print(message, highlight=False)


def main():
    """
    Main entrypoint for echosv subcommands
    """
    parser = ArgumentParser(prog="echosv", description=USAGE,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=["chain", "genotype", "match"], type=str, default=None,
                        help="Command to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit()
    args = parser.parse_args()
    if args.cmd == "chain":
        chain_main(" ".join(args.options))
    elif args.cmd == "genotype":
        genotype_main(" ".join(args.options))
    elif args.cmd == "match":    
        match_main(" ".join(args.options))


if __name__ == '__main__':
    main()