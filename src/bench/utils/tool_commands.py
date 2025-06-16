import json
import argparse
import sys
from pathlib import Path
from .functions import get_parse_function

def populate_tools(args, spacer_lendf=None):
    """
    Populate tools dictionary by reading all jsons from the tool_configs dir.
    """
    threads = args.threads
    contigs_file = f"{args.results_dir}/simulated_data/simulated_contigs.fa"
    spacers_file = f"{args.results_dir}/simulated_data/simulated_spacers.fa"
    results_dir = args.results_dir
    output_dir = f"{args.results_dir}/raw_outputs/"
    
    # Create the tool_configs directory if it doesn't exist
    config_dir = Path(__file__).parent.parent.parent.parent / "tool_configs"
    config_dir.mkdir(exist_ok=True)
    
    # Check if there are tool configurations available
    tool_configs = list(config_dir.glob("*.json"))
    # Load tool configurations from JSON files
    tools = {}
    for config_file in tool_configs:
        try:
            with open(config_file, 'r') as f:
                tool_config = json.load(f)
            
            # Make sure required keys are present
            required_keys = ['name', 'output_file', 'parse_function_name', 'script_name', 'command']
            if all(key in tool_config for key in required_keys):
                # Replace placeholders in the configuration
                tool_config['output_file'] = tool_config['output_file'].format(
                    output_dir=output_dir,
                    results_dir=results_dir
                )
                
                # Handle command differently based on whether it's a string or list
                if isinstance(tool_config['command'], list):
                    # Replace placeholders in command list or string
                    tool_config['command'] = [
                        cmd.format(
                            threads=threads,
                            contigs_file=contigs_file,
                            spacers_file=spacers_file,
                            output_dir=output_dir,
                            results_dir=results_dir
                        ) for cmd in tool_config['command']
                    ]
                
                # Add the tool to the dictionary
                tools[tool_config['name']] = tool_config
            else:
                print(f"Warning: Skipping tool config {config_file} due to missing required keys")
        
        except json.JSONDecodeError:
            print(f"Warning: Failed to load tool config from {config_file} - invalid JSON")
        except Exception as e:
            print(f"Warning: Error loading tool config from {config_file}: {str(e)}")
    
    # Adjust parse function based on spacer length dataframe
    for tool in tools.values():
        # If parse_function_name is not specified, default to parse_sam
        if 'parse_function_name' not in tool or not tool['parse_function_name']:
            tool['parse_function_name'] = 'parse_sam'
                
        # Add the actual function reference
        tool['parse_function'] = get_parse_function(tool['parse_function_name'])
    
    return tools

def prompt_for_tool_info():
    """Interactively prompt for tool information"""
    tool_info = {}
    
    print("Enter information for the new tool:")
    tool_info['name'] = input("Tool name (must be unique): ").strip()
    
    # Define output file
    default_output = f"{{output_dir}}/{tool_info['name']}_output.tsv"
    print(f"Default output file: {default_output}")
    custom_output = input(f"Custom output file (press Enter for default): ").strip()
    tool_info['output_file'] = custom_output if custom_output else default_output
    
    # Determine file type based on extension
    file_ext = tool_info['output_file'].split('.')[-1].lower()
    
    # Suggest appropriate parse function based on output type
    suggested_parse = "parse_tab"
    if file_ext == "sam":
        suggested_parse = "parse_samVn_with_lens_pysam"
    elif file_ext == "tsv" or file_ext == "blast6":
        suggested_parse = "parse_blastn_custom"
    
    print(f"Suggested parse function: {suggested_parse}")
    parse_func = input(f"Parse function name (press Enter for suggested): ").strip()
    tool_info['parse_function_name'] = parse_func if parse_func else suggested_parse
    
    # Script name
    default_script = f"{tool_info['name']}.sh"
    print(f"Default script name: {default_script}")
    custom_script = input(f"Custom script name (press Enter for default): ").strip()
    tool_info['script_name'] = custom_script if custom_script else default_script
    
    # Mamba environment
    tool_info['mamba_env'] = input(f"Mamba environment (leave empty for None): ").strip() or None
    
    # Command
    print("Now enter the command to run the tool.")
    print("You can use the following placeholders:")
    print("  {threads} - Number of threads")
    print("  {contigs_file} - Path to contigs file")
    print("  {spacers_file} - Path to spacers file")
    print("  {output_dir} - Path to output directory")
    print("  {results_dir} - Path to results directory")
    
    commands = []
    print("Enter command (one per line, empty line to finish):")
    while True:
        cmd = input("> ").strip()
        if cmd:
            commands.append(cmd.replace("{threads}", "str(threads)")
                               .replace("{contigs_file}", "contigs_file")
                               .replace("{spacers_file}", "spacers_file")
                               .replace("{output_dir}", "output_dir")
                               .replace("{results_dir}", "results_dir"))
        else:
            break
    
    if len(commands) == 1:
        tool_info['command'] = commands
    elif len(commands) > 1:
        # Join multiple commands with '&&'
        joined_cmd = " && ".join(commands)
        tool_info['command'] = [joined_cmd]
    else:
        print("No command entered. Using a placeholder command.")
        tool_info['command'] = ["echo 'No command specified'"]
    
    return tool_info

def add_tool_to_configs(tool_info):
    """Add a new tool configuration to the tool_configs directory"""
    # Get the directory where tool configs are stored
    config_dir = Path(__file__).parent.parent / "tool_configs"
    config_dir.mkdir(exist_ok=True)
    
    # Check if a tool with this name already exists
    config_path = config_dir / f"{tool_info['name']}.json"
    if config_path.exists():
        overwrite = input(f"A tool with the name '{tool_info['name']}' already exists. Overwrite? (y/n): ").strip().lower()
        if overwrite != 'y':
            print("Tool addition canceled.")
            return False
    
    # Write the tool configuration to a file
    with open(config_path, 'w') as f:
        json.dump(tool_info, f, indent=2)
    
    print(f"Successfully added tool '{tool_info['name']}' to {config_path}")
    return True

def list_available_tools():
    """List all available tools from the tool_configs directory"""
    config_dir = Path(__file__).parent.parent / "tool_configs"
    if not config_dir.exists():
        print("No tool configurations available.")
        return
    
    tools = list(config_dir.glob("*.json"))
    if not tools:
        print("No tool configurations available.")
        return
    
    print(f"Available tools ({len(tools)}):")
    for tool_path in tools:
        try:
            with open(tool_path, 'r') as f:
                tool_config = json.load(f)
            print(f" - {tool_config['name']} ({tool_path.name})")
        except:
            print(f" - Unknown tool ({tool_path.name}) - invalid configuration")

def main():
    parser = argparse.ArgumentParser(description="Tool commands management")
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Add tool subcommand
    add_parser = subparsers.add_parser("add-tool", help="Add a new tool")
    add_parser.add_argument("--mode", choices=["interactive", "json"], default="interactive", 
                          help="Mode to add a tool: interactive prompts or from JSON file")
    add_parser.add_argument("--json-file", help="Path to JSON file with tool configuration")
    
    # List tools subcommand
    list_parser = subparsers.add_parser("list-tools", help="List all available tools")
    
    # Convert existing tools subcommand
    convert_parser = subparsers.add_parser("convert-default-tools", 
                                         help="Create JSON files for default tools")
    
    args = parser.parse_args()
    
    if args.command == "add-tool":
        if args.mode == "interactive":
            tool_info = prompt_for_tool_info()
            if tool_info:
                add_tool_to_configs(tool_info)
        
        elif args.mode == "json":
            if not args.json_file:
                print("Error: --json-file is required with --mode json")
                sys.exit(1)
            
            try:
                with open(args.json_file, 'r') as f:
                    tool_info = json.load(f)
                
                required_keys = ['name', 'output_file', 'parse_function_name', 'script_name', 'command']
                missing_keys = [key for key in required_keys if key not in tool_info]
                
                if missing_keys:
                    print(f"Error: JSON file is missing required keys: {', '.join(missing_keys)}")
                    sys.exit(1)
                
                add_tool_to_configs(tool_info)
            
            except json.JSONDecodeError:
                print(f"Error: {args.json_file} is not a valid JSON file")
                sys.exit(1)
            except FileNotFoundError:
                print(f"Error: File {args.json_file} not found")
                sys.exit(1)
    
    elif args.command == "list-tools":
        list_available_tools()
        
    else:
        parser.print_help()

if __name__ == "__main__":
    main() 