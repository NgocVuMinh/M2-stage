# Analyze multi-site sampling within dt.ebi2.p
multi_site_analysis <- sel %>%
  # Group by patient, tissue, and stage
  group_by(Sample.Characteristic.individual., 
           Sample.Characteristic.organism.part., pcw_stage
           ) %>%
  # Summarize sampling sites for each combination
  summarise(
    n_samples = n(),
    n_sites = n_distinct(Sample.Characteristic.sampling.site.),
    sampling_sites = paste(unique(Sample.Characteristic.sampling.site.), collapse = ", "),
    run_ids = paste(Run, collapse = ", "),
    .groups = 'drop'
  ) %>%
  # Filter for cases with multiple sites
  filter(n_sites > 1)

# Display patients with multi-site sampling
print(multi_site_analysis)

# Save to dataframe
patients_multisite <- multi_site_analysis

# Analyze which tissues are commonly sampled at multiple sites
tissue_multisite_summary <- multi_site_analysis %>%
  group_by(Sample.Characteristic.organism.part.) %>%
  summarise(
    n_patients_with_multisite = n_distinct(Sample.Characteristic.individual.),
    n_total_multisite_instances = n(),
    stages_involved = paste(unique(pcw_stage), collapse = ", "),
    site_combinations = paste(unique(sampling_sites), collapse = " | "),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_total_multisite_instances))

print("\nTissues commonly sampled at multiple sites:")
print(tissue_multisite_summary)

# Analyze which stages have multi-site sampling
stage_multisite_summary <- multi_site_analysis %>%
  group_by(pcw_stage) %>%
  summarise(
    n_patients_with_multisite = n_distinct(Sample.Characteristic.individual.),
    n_total_multisite_instances = n(),
    tissues_involved = paste(unique(Sample.Characteristic.organism.part.), collapse = ", "),
    site_combinations = paste(unique(sampling_sites), collapse = " | "),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_total_multisite_instances))

print("\nStages with multi-site sampling:")
print(stage_multisite_summary)

# Analyze tissue-stage combinations with multi-site sampling
tissue_stage_multisite <- multi_site_analysis %>%
  group_by(Sample.Characteristic.organism.part., pcw_stage) %>%
  summarise(
    n_patients = n_distinct(Sample.Characteristic.individual.),
    n_instances = n(),
    site_combinations = paste(unique(sampling_sites), collapse = " | "),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_instances))

print("\nTissue-Stage combinations with multi-site sampling:")
print(tissue_stage_multisite)