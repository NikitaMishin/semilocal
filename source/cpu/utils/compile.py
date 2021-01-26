def construct_compile_commands(folder, compiler_path):
    """
    return commands to run in shell to prepare implementations to run
    :return:
    """
    rm_and_mkdirs = [
        ['rm', '-rf', folder],
        ['mkdir', folder]
    ]

    cmake_command = ['cmake', f'-B{folder}',
                     '-H.', '-DCMAKE_BUILD_TYPE=Release',
                     f'-DCMAKE_CXX_COMPILER={compiler_path}']

    make = ['make', '-C', folder]

    rm_and_mkdirs.extend([cmake_command, make])

    return rm_and_mkdirs
