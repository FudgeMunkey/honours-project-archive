# Go through the GSW documentation and pull out the functions
# Stub them for easy hypothesis testing

DOCUMENTATION_PATH = "all_functions.txt"


def extract_all_functions():
    all_parameters = []
    full_functions = {}

    f = open(DOCUMENTATION_PATH, "r")
    lines = f.read().split("\n")

    for line in lines:
        if "gsw." in line:
            # Extract the function name
            function_name = line[4 : line.index("(")]

            # Extract the parameters
            raw_parameters = line[line.index("(") + 1 : line.index(")")].split(", ")
            parameters = []

            for parameter in raw_parameters:
                if "=" in parameter:
                    parameter = parameter[: parameter.index("=")]

                parameters.append(parameter)

            # Track new strategies
            for parameter in parameters:
                new_parameter = parameter + "_st"

                if new_parameter not in all_parameters:
                    all_parameters.append(new_parameter)

            # Build the function
            given_parameters = ["\t" + parameter + "_st" for parameter in parameters]
            given_parameters = ",\n".join(given_parameters) + ","
            given_text = f"@given(\n{given_parameters}\n)"

            function_parameters = "self, " + ", ".join(parameters)
            function_declaration = f"def test_{function_name}({function_parameters}):\n\tresult = {function_name}({function_parameters})\n\n"

            full_function = f"{given_text}\n{function_declaration}"
            full_functions[function_name] = full_function

    return all_parameters, full_functions


def extract_subpackage_function_names():
    function_subpackage = {}

    subpackages = [
        "ConversionFunctions",
        "Density",
        "Energy",
        "Stability",
        "Geostrophy",
        "Ice",
        "Freezing",
    ]

    for subpackage in subpackages:
        f = open(f"subpackages/{subpackage}.txt")
        lines = f.read().split("\n")

        for line in lines:
            if "gsw." in line:
                function_name = line[line.index(".", 4) + 1 : line.index("(")]
                function_subpackage[function_name] = subpackage

    return function_subpackage


if __name__ == "__main__":
    function_subpackage = extract_subpackage_function_names()
    all_parameters, full_functions = extract_all_functions()

    function_classes = {
        "ConversionFunctions": [],
        "Density": [],
        "Energy": [],
        "Stability": [],
        "Geostrophy": [],
        "Ice": [],
        "Freezing": [],
        "Misc": [],
    }

    # Print Results
    print("# STRATEGIES")
    for parameter in all_parameters:
        print(f"{parameter} = None")

    print("\n")

    for function in full_functions:
        full_function = full_functions[function]

        if function in function_subpackage:
            subpackage = function_subpackage[function]
            function_classes[subpackage].append(full_function)
        else:
            function_classes["Misc"].append(full_function)

    for function_class in function_classes:
        print(f"class Test{function_class}:")

        full_class_functions = function_classes[function_class]

        for full_function in full_class_functions:
            full_function = "\n".join(["\t" + ff for ff in full_function.split("\n")])
            print(full_function)
