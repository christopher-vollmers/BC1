use pyo3::PyResult;
use pyo3::types::PyDict;
use std::collections::HashMap;
use std::fs::File;
use needletail::{parse_fastx_file, FastxReader};
use pyo3::prelude::*;
use pyo3::{pymodule, wrap_pyfunction, pyfunction, Python};
use pyo3::types::PyModule;

use std::io::prelude::*;

use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;


#[pyfunction]
fn create_labeled_subreads(UMIdict: HashMap<String, (String, u32)>, subread_files: &str, out2: &str) -> PyResult<()> {
    let mut out_sub = File::create(out2).expect("Failed to create output file");
    
    let files: Vec<&str> = subread_files.split(',').collect();
    
    for subread_file in files {
        println!("Parsing file start!");
        let mut reader = parse_fastx_file(subread_file).expect("Expecting a valid FASTX file");
        println!("File parsed!");

        while let Some(record) = reader.next() {
            let seq_record = record.expect("Failed to read record");
            let name = std::str::from_utf8(seq_record.id().as_ref()).expect("Failed to convert sequence ID").to_string();
            let seq = std::str::from_utf8(seq_record.seq().as_ref()).expect("Failed to convert sequence").to_string();
            let q = std::str::from_utf8(seq_record.qual().unwrap().as_ref()).expect("Failed to convert quality scores").to_string();
            
            let root = name.split('_').next().unwrap();
            
            if let Some((UMI, number)) = UMIdict.get(root) {
                let line = format!("{}\t{}\t{}\t{}\t{}\n", number, UMI, name, seq, q);
                out_sub.write_all(line.as_bytes()).expect("Failed to write to output file");
            }
        }
    }
    
    out_sub.flush().expect("Failed to flush output file");
    Ok(())
}

#[pyfunction]
fn extracting_umis_starsolo(py:Python, input_sam: &str, output_file_root: &str, umi_patterns: &str) -> PyResult<(PyObject, i32)> {
    println!("Starting to extract UMIs in rust");
    let input_file = File::open(input_sam)?;
    let reader = BufReader::new(input_file);
    let mut output_file = File::create(format!("{}.UMIs", output_file_root))?;

    let mut umi_dict: HashMap<String, (String, usize)> = HashMap::new();
    let median_umi_length: f32 = 0.0; 
    let mut umi_numbers = HashMap::new();
    let mut counter = 0;
    let mut reason_dict = HashMap::new();
    let mut umi_lengths = Vec::new();

    let umi_patterns: Vec<&str> = umi_patterns.split(',').collect();


    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        let name = parts[0];
        let tags = &parts[8..];
        let mut umis = Vec::new();
    
        for &umi_pattern in &umi_patterns {
            let mut umi = tags.iter()
                .find(|&&tag| tag.starts_with(umi_pattern))
                .map_or(String::new(), |&tag| tag.replacen(umi_pattern, "", 1));
    
            if umi.ends_with("-1") {
                umi.pop();
                umi.pop();
            }
            umi = umi.replace('_', "");
    
            let reason = if umi.is_empty() {
                format!("no {} found", umi_pattern)
            } else {
                format!("{} found", umi_pattern)
            };
            *reason_dict.entry(reason).or_insert(0) += 1;
            umis.push(umi);
        }

        let umi = umis.join("_");
        if !umi_numbers.contains_key(&umi) {
            counter += 1;
            umi_numbers.insert(umi.clone(), counter);
        }
        writeln!(output_file, "{}\t{}", name, umi)?;

        let root = name.split('_').next().unwrap_or_default().to_string();
        umi_dict.insert(root.clone(), (umi.clone(), *umi_numbers.get(&umi).unwrap()));
        umi_lengths.push(umi.len() as i32);
    }

    let median_umi_length = median(&umi_lengths);

    let umi_dict_py = PyDict::new(py);
    for (key, (umi, count)) in umi_dict.iter() {
        umi_dict_py.set_item(key, (umi, count))?;
    }

    let median_umi_length_int: i32 = median_umi_length.trunc() as i32;

    Ok((umi_dict_py.to_object(py), median_umi_length_int))
}

fn median(numbers: &[i32]) -> f32 {
    let mut numbers = numbers.to_vec();
    numbers.sort_unstable();
    let mid = numbers.len() / 2;
    if numbers.len() % 2 == 0 {
        (numbers[mid - 1] as f32 + numbers[mid] as f32) / 2.0
    } else {
        numbers[mid] as f32
    }
}

#[pymodule]
fn rust_tools(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(create_labeled_subreads, m)?)?;
    m.add_function(wrap_pyfunction!(extracting_umis_starsolo, m)?)?;
    Ok(())
}
