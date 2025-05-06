import subprocess
import tempfile
import polars as pl

def parse_time_to_seconds(time_str):
    """Convert Slurm time format to seconds"""
    if not time_str or time_str == "":
        return 0
    try:
        # Handle days-hours:minutes:seconds format
        if '-' in time_str:
            days, time_part = time_str.split('-')
            hours, minutes, seconds = time_part.split(':')
            return (int(days) * 24 * 3600 + 
                   int(hours) * 3600 + 
                   int(minutes) * 60 + 
                   int(seconds))
        # Handle hours:minutes:seconds format
        else:
            parts = time_str.split(':')
            if len(parts) == 3:
                hours, minutes, seconds = parts
                return int(hours) * 3600 + int(minutes) * 60 + int(seconds)
            return 0
    except (ValueError, AttributeError):
        return 0

def parse_mem(mem_str):
    if not mem_str:
        return 0
    try:
        if mem_str.endswith('K'):
            return float(mem_str[:-1]) / 1024
        elif mem_str.endswith('G'):
            return float(mem_str[:-1]) * 1024
        elif mem_str.endswith('T'):
            return float(mem_str[:-1]) * 1024 * 1024
        else:
            return float(mem_str)
    except (ValueError, AttributeError):
        return 0

def format_seconds(seconds):
    days = seconds // (24 * 3600)
    remaining = seconds % (24 * 3600)
    hours = remaining // 3600
    remaining %= 3600
    minutes = remaining // 60
    seconds = remaining % 60
    if days > 0:
        return f"{days}-{hours:02d}:{minutes:02d}:{seconds:02d}"
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"

def pyseff(sacct_df=None):
    """
    Main function to analyze job efficiency from sacct output, similar to Slurm's seff
    
    Uses specific fields from sacct:
    - JobID: for identifying jobs
    - JobName: name of the job
    - State: job state
    - ExitCode: exit status
    - AllocCPUS: allocated CPUs
    - Elapsed: total job time
    - TotalCPU: total CPU time (user + system)
    - MaxRSS: maximum resident set size
    - ReqMem: requested memory
    """
    if sacct_df is None:
        # Create temporary file for sacct output with specific fields
        tmp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
        sacct_cmd = (
            "sacct -u $USER --parsable "
            "--format=JobID,JobName,State,ExitCode,AllocCPUS,Elapsed,"
            "TotalCPU,MaxRSS,ReqMem "
            "--starttime 2025-01-01 "
            "> " + tmp_file.name
        )
        subprocess.run(sacct_cmd, shell=True, text=True, capture_output=True)
        sacct_df = pl.read_csv(tmp_file.name, separator='|')
        tmp_file.close()
    sacct_df = sacct_df.with_columns(
        pl.col('JobID').str.split('.').list.first().alias('BaseJobID')
    )
    cancelled_jobs = sacct_df.filter(pl.col('State').str.contains_any(patterns=['CANCELLED']))['BaseJobID'].unique().to_list()
    
    # Get base job IDs and filter out cancelled jobs
    sacct_df = sacct_df.filter(~pl.col('BaseJobID').is_in(cancelled_jobs))
    
    # Calculate CPU times in seconds
    sacct_df = sacct_df.with_columns([
        pl.col('TotalCPU').map_elements(parse_time_to_seconds,return_dtype=pl.UInt32).alias('TotalCPU_Seconds'),
        pl.col('Elapsed').map_elements(parse_time_to_seconds,return_dtype=pl.UInt32).alias('Elapsed_Seconds')
    ])
    
    # Aggregate by BaseJobID
    sacct_df = sacct_df.group_by('BaseJobID').agg([
        pl.col('JobName').first().alias('JobName'),
        pl.col('AllocCPUS').max().alias('AllocCPUS'),
        pl.col('State').unique().alias('State'),
        pl.col('ExitCode').unique().alias('ExitCode'),
        pl.col('MaxRSS').max().alias('MaxRSS'),
        pl.col('ReqMem').first().alias('ReqMem'),
        pl.col('Elapsed_Seconds').sum().alias('Elapsed_Seconds'),
        pl.col('TotalCPU_Seconds').sum().alias('TotalCPU_Seconds')
    ])
    
    # Calculate CPU efficiency
    sacct_df = sacct_df.with_columns([
        (pl.col('TotalCPU_Seconds') / 
         (pl.col('Elapsed_Seconds') * pl.col('AllocCPUS')) * 100)
        .round(2)
        .alias('CPU_Efficiency')
    ])
    
    # Calculate memory efficiency
    sacct_df = sacct_df.with_columns([
        pl.col('MaxRSS').map_elements(parse_mem,return_dtype=pl.Float64).alias('MaxRSS_MB'),
        pl.col('ReqMem').map_elements(parse_mem,return_dtype=pl.Float64).alias('ReqMem_MB'),
        
        (pl.col('MaxRSS').map_elements(parse_mem,return_dtype=pl.Float64) / 
         pl.col('ReqMem').map_elements(parse_mem,return_dtype=pl.Float64) * 100)
        .round(2)
        .alias('Memory_Efficiency')
    ])
    
    # Format time columns back to human-readable format
    result_df = sacct_df.with_columns([
        pl.col('Elapsed_Seconds').map_elements(format_seconds,return_dtype=pl.Utf8).alias('Elapsed'),
        pl.col('TotalCPU_Seconds').map_elements(format_seconds,return_dtype=pl.Utf8).alias('TotalCPU')
    ]).select([
        'BaseJobID',
        'JobName',
        'State',
        'ExitCode',
        'AllocCPUS',
        'Elapsed',
        'TotalCPU',
        'CPU_Efficiency',
        'MaxRSS_MB',
        'ReqMem_MB',
        'Memory_Efficiency'
    ])
    
    return result_df
