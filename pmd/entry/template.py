import inquirer


def main():
    questions = [
        inquirer.List(
            'system',
            message="What system do you need?",
            choices=[
                'System (amorphous homopolymer)',
                'SolventSystem (homopolymer + solvent)',
                'GasSystem (homopolymer + gas)'
            ],
        ),
        inquirer.List(
            'builder',
            message=
            "What force field (builder) do you want to use for this system?",
            choices=[
                'opls-lbcc (PSP)',
                'opls-cm1a (PSP)',
                'gaff2-gasteiger (PSP)',
                'gaff2-am1bcc (PSP)',
                'pcff (EMC)',
                'opls-aa (EMC)',
                'opls-ua (EMC)',
                'trappe (EMC)',
            ],
        ),
        inquirer.List(
            'lammps',
            message="What property do you want to compute?",
            choices=[
                'Glass transition temperature', 'Gas/solvent diffusivity',
                'Viscosity', 'Mechanical properties', 'Thermal conductivity'
            ],
        ),
        inquirer.List(
            'job',
            message="What job scheduler system do you use?",
            choices=['Torque', 'Slurm', 'N/A (run locally)'],
        ),
    ]
    answers = inquirer.prompt(questions)
    print(answers)


if __name__ == '__main__':
    main()
