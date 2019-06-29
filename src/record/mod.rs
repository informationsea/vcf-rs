use super::*;

/// A VCF record structure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VCFRecord {
    pub chromosome: String,
    pub position: u64,
    pub id: Vec<String>,
    pub reference: String,
    pub alternative: Vec<String>,
    pub quality: Option<String>,
    pub filter: Vec<String>,
    pub info: IndexMap<String, Vec<String>>,
    pub format: Vec<String>,
    pub call: HashMap<String, HashMap<String, Vec<String>>>,
}

impl VCFRecord {
    /// Create VCFRecord from a text line.    
    pub fn parse_line(line: &str, samples: &[String]) -> Result<VCFRecord, VCFParseError> {
        let elements: Vec<_> = line.trim().split('\t').collect();

        if elements.len() < 5 {
            return Err(VCFParseError::NotEnoughColumns);
        }

        let info = if elements.len() < 8 || elements[7] == "." {
            IndexMap::new()
        } else {
            elements[7]
                .split(';')
                .map(|x| {
                    let mut y = x.splitn(2, '=');
                    let first = y.next().unwrap();
                    let second = y.next();
                    (
                        first.to_string(),
                        second
                            .map(|x| {
                                x.split(',')
                                    .map(|z| decode_vcf_value(z))
                                    .collect::<Vec<_>>()
                            })
                            .unwrap_or_else(|| vec![]),
                    )
                })
                .collect()
        };

        let format = if elements.len() < 9 || elements[8] == "." {
            vec![]
        } else {
            elements[8]
                .split(':')
                .map(|x| decode_vcf_value(x))
                .collect()
        };

        let call = elements
            .iter()
            .skip(9)
            .zip(samples)
            .map(|(v, k)| {
                (
                    k.to_string(),
                    v.split(':')
                        .zip(format.iter())
                        .filter(|(x, _)| *x != ".")
                        .map(|(x, f)| {
                            (
                                f.to_string(),
                                x.split(',').map(|z| decode_vcf_value(z)).collect(),
                            )
                        })
                        .collect(),
                )
            })
            .collect();

        Ok(VCFRecord {
            chromosome: decode_vcf_value(elements[0]),
            position: elements[1]
                .parse::<u64>()
                .map_err(|_| VCFParseError::PositionIsNotNumber(elements[1].to_string()))?,
            id: if elements[2] == "." {
                Vec::new()
            } else {
                elements[2]
                    .split(';')
                    .map(|x| decode_vcf_value(x))
                    .collect()
            },
            reference: decode_vcf_value(elements[3]),
            alternative: elements[4]
                .split(',')
                .map(|x| decode_vcf_value(x))
                .collect(),
            quality: dot_value(&elements, 5),
            filter: if elements.len() < 7 || elements[6] == "." {
                vec![]
            } else {
                elements[6]
                    .split(',')
                    .map(|x| decode_vcf_value(x))
                    .collect()
            },
            info,
            format,
            call,
        })
    }

    /// write a VCF line.
    pub fn write_line<W: io::Write>(&self, writer: &mut W, samples: &[String]) -> io::Result<()> {
        write!(
            writer,
            "{}\t{}\t",
            encode_vcf_value(&self.chromosome),
            self.position,
        )?;
        if self.id.is_empty() {
            writer.write_all(b".")?;
        } else {
            for (i, one) in self.id.iter().enumerate() {
                if i != 0 {
                    writer.write_all(b";")?;
                }
                writer.write_all(encode_vcf_value(one).as_bytes())?;
            }
        }
        write!(writer, "\t{}\t", self.reference)?;
        write_encoded_vec_or_dot(writer, &self.alternative)?;
        writer.write_all(b"\t")?;
        write_str_or_dot(writer, &self.quality)?;
        writer.write_all(b"\t")?;
        write_encoded_vec_or_dot(writer, &self.filter)?;

        writer.write_all(b"\t")?;
        if self.info.is_empty() {
            writer.write_all(b".")?;
        } else {
            for (i, (k, v)) in self.info.iter().enumerate() {
                if i != 0 {
                    writer.write_all(b";")?;
                }
                writer.write_all(k.as_bytes())?;
                if v.is_empty() {
                    // skip
                } else {
                    writer.write_all(b"=")?;
                    write_encoded_vec_or_dot(writer, &v)?;
                }
            }
        }

        if !samples.is_empty() {
            writer.write_all(b"\t")?;
            if self.format.is_empty() {
                writer.write_all(b".")?;
            } else {
                for (i, one) in self.format.iter().enumerate() {
                    if i != 0 {
                        writer.write_all(b":")?;
                    }
                    writer.write_all(encode_vcf_value(&one).as_bytes())?;
                }
            }
            writer.write_all(b"\t")?;
            for (si, one_sample) in samples.iter().enumerate() {
                if si != 0 {
                    writer.write_all(b"\t")?;
                }
                if let Some(call_result) = self.call.get(one_sample) {
                    for (i, one_format) in self.format.iter().enumerate() {
                        if i != 0 {
                            writer.write_all(b":")?;
                        }
                        if let Some(one_call) = call_result.get(one_format) {
                            write_encoded_vec_or_dot(writer, one_call)?;
                        } else {
                            writer.write_all(b".")?;
                        }
                    }
                } else {
                    writer.write_all(b".")?;
                }
            }
        }

        writer.write_all(b"\n")?;

        Ok(())
    }
}

#[cfg(test)]
mod test;