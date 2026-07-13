//! `sam-diff <oracle.sam> <ours.sam> [--json <path>]`
//!
//! Prints a JSON concordance report (and optionally writes it to a file).

use std::path::Path;
use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let mut positional = Vec::new();
    let mut json_out: Option<String> = None;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--json" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--json requires a path");
                    return ExitCode::from(2);
                }
                json_out = Some(args[i].clone());
            }
            other => positional.push(other.to_string()),
        }
        i += 1;
    }
    if positional.len() != 2 {
        eprintln!("usage: sam-diff <oracle.sam> <ours.sam> [--json <path>]");
        return ExitCode::from(2);
    }

    let oracle = match bwa_diff::parse_primary(Path::new(&positional[0])) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("error reading {}: {e}", positional[0]);
            return ExitCode::from(2);
        }
    };
    let ours = match bwa_diff::parse_primary(Path::new(&positional[1])) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("error reading {}: {e}", positional[1]);
            return ExitCode::from(2);
        }
    };

    let report = bwa_diff::compare(&oracle, &ours, 20);
    let json = serde_json::to_string_pretty(&report).expect("serialize report");
    if let Some(path) = &json_out {
        if let Err(e) = std::fs::write(path, &json) {
            eprintln!("error writing {path}: {e}");
            return ExitCode::from(2);
        }
    }
    println!("{json}");
    ExitCode::SUCCESS
}
