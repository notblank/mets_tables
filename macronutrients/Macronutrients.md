```R
library(tidyverse)
library(broom)
library(car)
library(rms)
```

    ── [1mAttaching packages[22m ─────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──
    
    [32m✔[39m [34mggplot2[39m 3.3.3     [32m✔[39m [34mpurrr  [39m 0.3.4
    [32m✔[39m [34mtibble [39m 3.1.0     [32m✔[39m [34mdplyr  [39m 1.0.5
    [32m✔[39m [34mtidyr  [39m 1.1.3     [32m✔[39m [34mstringr[39m 1.4.0
    [32m✔[39m [34mreadr  [39m 1.4.0     [32m✔[39m [34mforcats[39m 0.5.1
    
    ── [1mConflicts[22m ────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m    masks [34mstats[39m::lag()
    
    Loading required package: carData
    
    
    Attaching package: ‘car’
    
    
    The following object is masked from ‘package:dplyr’:
    
        recode
    
    
    The following object is masked from ‘package:purrr’:
    
        some
    
    
    Loading required package: Hmisc
    
    Loading required package: lattice
    
    Loading required package: survival
    
    Loading required package: Formula
    
    
    Attaching package: ‘Hmisc’
    
    
    The following objects are masked from ‘package:dplyr’:
    
        src, summarize
    
    
    The following objects are masked from ‘package:base’:
    
        format.pval, units
    
    
    Loading required package: SparseM
    
    
    Attaching package: ‘SparseM’
    
    
    The following object is masked from ‘package:base’:
    
        backsolve
    
    
    
    Attaching package: ‘rms’
    
    
    The following objects are masked from ‘package:car’:
    
        Predict, vif
    
    



```R
macro_nut_portions_info.csv <- read_csv('../data/macro_nut_portions_info.csv')
conditions_info <- read_csv('../data/conditions_info.csv') 
```

    
    [36m──[39m [1m[1mColumn specification[1m[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    cols(
      id = [32mcol_double()[39m,
      carbohydrates = [32mcol_double()[39m,
      energy = [32mcol_double()[39m,
      fiber = [32mcol_double()[39m,
      lipids = [32mcol_double()[39m,
      protein = [32mcol_double()[39m,
      energy_macro = [32mcol_double()[39m,
      percent_carbohydrate = [32mcol_double()[39m,
      percent_protein = [32mcol_double()[39m,
      percent_lipids = [32mcol_double()[39m,
      beef_cat = [31mcol_character()[39m,
      dairy_cat = [31mcol_character()[39m,
      legumes_cat = [31mcol_character()[39m,
      white_meat_cat = [31mcol_character()[39m,
      beef = [32mcol_double()[39m,
      dairy = [32mcol_double()[39m,
      legumes = [32mcol_double()[39m,
      white_meat = [32mcol_double()[39m
    )
    
    
    
    [36m──[39m [1m[1mColumn specification[1m[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    cols(
      .default = col_double(),
      education = [31mcol_character()[39m,
      sex = [31mcol_character()[39m,
      location_type = [31mcol_character()[39m,
      smokes = [31mcol_character()[39m,
      drinks = [31mcol_character()[39m,
      elevated_bp = [33mcol_logical()[39m,
      elevated_triglycerids_mg_d_l = [33mcol_logical()[39m,
      elevated_glucose_mg_d_l = [33mcol_logical()[39m,
      low_hdl_mg_d_l = [33mcol_logical()[39m,
      met_s = [33mcol_logical()[39m,
      d_diagnosed = [33mcol_logical()[39m,
      h_diagnosed = [33mcol_logical()[39m
    )
    [36mℹ[39m Use [30m[47m[30m[47m`spec()`[47m[30m[49m[39m for the full column specifications.
    
    



```R
macro_info <- conditions_info %>% 
    inner_join(macro_nut_portions_info.csv, by = 'id')
```

## Macronutrients and conditions


```R
macro_info %>% 
    select(contains("percent_")) %>% 
    pivot_longer(cols = contains("percent_"), names_to = "source", values_to = "percent") %>% 
    group_by(source) %>% 
    summarise(m_macro = round(100*mean(percent), 1), sd_macro = round(100*sd(percent), 1))
```


<table class="dataframe">
<caption>A tibble: 3 × 3</caption>
<thead>
	<tr><th scope=col>source</th><th scope=col>m_macro</th><th scope=col>sd_macro</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>percent_carbohydrate</td><td>54.7</td><td>7.8</td></tr>
	<tr><td>percent_lipids      </td><td>22.8</td><td>5.3</td></tr>
	<tr><td>percent_protein     </td><td>22.5</td><td>4.9</td></tr>
</tbody>
</table>




```R
t_test_macro <- function(condition, macro){
    
    t_df <- macro_info %>% 
                mutate(overweight_obesity = ifelse(bmi >= 25, TRUE, FALSE)) %>% 
                na.omit() %>% 
                select(condition, macro) 

    t <- t.test(x = t_df %>% filter(get(condition) == TRUE), 
                y = t_df %>% filter(get(condition) == FALSE))
    
    return(t)
}
```


```R
energy_met_s <- macro_info %>% 
    na.omit() %>% 
    select(met_s, energy) %>% 
    pivot_longer(cols = -met_s, names_to = 'source', values_to = 'percent') %>% 
    group_by(met_s, source) %>% 
    summarise(m_percent = mean(percent), sd_percent = sd(percent)) %>% 
    mutate(m_percent = round(m_percent, 0), sd_percent = round(sd_percent, 0), 
           dist = paste(m_percent,sd_percent, sep='')) %>% 
    select(met_s, source, dist) %>% 
    pivot_wider(names_from = met_s, values_from = dist) %>% 
    mutate(condition = 'met_s') %>% 
    bind_cols(p_val = round(t_test_macro("met_s", "energy")$p.val, 3))
```

    `summarise()` has grouped output by 'met_s'. You can override using the `.groups` argument.
    
    Note: Using an external vector in selections is ambiguous.
    [34mℹ[39m Use `all_of(condition)` instead of `condition` to silence this message.
    [34mℹ[39m See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    [90mThis message is displayed once per session.[39m
    
    Note: Using an external vector in selections is ambiguous.
    [34mℹ[39m Use `all_of(macro)` instead of `macro` to silence this message.
    [34mℹ[39m See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    [90mThis message is displayed once per session.[39m
    



```R
mets_macro <- macro_info %>% 
    na.omit() %>% 
    select(met_s, contains('percent')) %>% 
    pivot_longer(cols = -met_s, names_to = 'source', values_to = 'percent') %>% 
    mutate(source = str_remove(source, 'percent_')) %>% 
    group_by(met_s, source) %>% 
    summarise(m_percent = mean(percent), sd_percent = sd(percent)) %>% 
    mutate(m_percent = round(100*m_percent, 0), sd_percent = round(100*sd_percent, 0), 
           dist = paste(m_percent, '% ', sd_percent, '%', sep='')) %>% 
    select(met_s, source, dist) %>% 
    pivot_wider(names_from = met_s, values_from = dist) %>% 
    mutate(condition = 'met_s') %>% 
    bind_cols(tibble(p_val = c(t_test_macro("met_s", "carbohydrates")$p.val, 
                               t_test_macro("met_s", "lipids")$p.val, 
                               t_test_macro("met_s", "protein")$p.val))) %>% 
    mutate(p_val = round(p_val, 3))
```

    `summarise()` has grouped output by 'met_s'. You can override using the `.groups` argument.
    



```R
dm2_macro <- macro_info %>% 
    na.omit() %>% 
    select(d_diagnosed, contains('percent')) %>% 
    pivot_longer(cols = -d_diagnosed, names_to = 'source', values_to = 'percent') %>% 
    mutate(source = str_remove(source, 'percent_')) %>% 
    group_by(d_diagnosed, source) %>% 
    summarise(m_percent = mean(percent), sd_percent = sd(percent)) %>% 
    mutate(m_percent = round(100*m_percent, 0), sd_percent = round(100*sd_percent, 0), 
           dist = paste(m_percent, '% ', sd_percent, '%', sep='')) %>% 
    select(d_diagnosed, source, dist) %>% 
    pivot_wider(names_from = d_diagnosed, values_from = dist)%>% 
    mutate(condition = 'dm2') %>% 
    bind_cols(tibble(p_val = c(t_test_macro("d_diagnosed", "carbohydrates")$p.val, 
                               t_test_macro("d_diagnosed", "lipids")$p.val, 
                               t_test_macro("d_diagnosed", "protein")$p.val))) %>% 
    mutate(p_val = round(p_val, 3))
```

    `summarise()` has grouped output by 'd_diagnosed'. You can override using the `.groups` argument.
    



```R
energy_dm2 <- macro_info %>% 
    na.omit() %>% 
    select(d_diagnosed, energy) %>% 
    pivot_longer(cols = -d_diagnosed, names_to = 'source', values_to = 'percent') %>% 
    group_by(d_diagnosed, source) %>% 
    summarise(m_percent = mean(percent), sd_percent = sd(percent)) %>% 
    mutate(m_percent = round(m_percent, 0), sd_percent = round(sd_percent, 0), 
           dist = paste(m_percent, sd_percent, sep='')) %>% 
    select(d_diagnosed, source, dist) %>% 
    pivot_wider(names_from = d_diagnosed, values_from = dist) %>% 
    mutate(condition = 'd_diagnosed') %>% 
    bind_cols(p_val = round(t_test_macro("d_diagnosed", "energy")$p.val, 3))
```

    `summarise()` has grouped output by 'd_diagnosed'. You can override using the `.groups` argument.
    



```R
overob_macro <- macro_info %>% 
    mutate(overweight_obesity = ifelse(bmi >= 25, TRUE, FALSE)) %>% 
    na.omit() %>% 
    select(overweight_obesity, contains('percent')) %>% 
    pivot_longer(cols = -overweight_obesity, names_to = 'source', values_to = 'percent') %>% 
    mutate(source = str_remove(source, 'percent_')) %>% 
    group_by(overweight_obesity, source) %>% 
    summarise(m_percent = mean(percent), sd_percent = sd(percent)) %>% 
    mutate(m_percent = round(100*m_percent, 0), sd_percent = round(100*sd_percent, 0), 
           dist = paste(m_percent, '% ', sd_percent, '%', sep='')) %>% 
    select(overweight_obesity, source, dist) %>% 
    pivot_wider(names_from = overweight_obesity, values_from = dist) %>% 
    mutate(condition = 'overweight_obesity') %>% 
    bind_cols(tibble(p_val = c(t_test_macro("overweight_obesity", "carbohydrates")$p.val, 
                               t_test_macro("overweight_obesity", "lipids")$p.val, 
                               t_test_macro("overweight_obesity", "protein")$p.val))) %>% 
    mutate(p_val = round(p_val, 3))
```

    `summarise()` has grouped output by 'overweight_obesity'. You can override using the `.groups` argument.
    



```R
energy_overob <- macro_info %>% 
    mutate(overweight_obesity = ifelse(bmi >= 25, TRUE, FALSE)) %>% 
    na.omit() %>% 
    select(overweight_obesity, energy) %>% 
    pivot_longer(cols = -overweight_obesity, names_to = 'source', values_to = 'percent') %>% 
    group_by(overweight_obesity, source) %>% 
    summarise(m_percent = mean(percent), sd_percent = sd(percent)) %>% 
    mutate(m_percent = round(m_percent, 0), sd_percent = round(sd_percent, 0), 
           dist = paste(m_percent, sd_percent, sep='')) %>% 
    select(overweight_obesity, source, dist) %>% 
    pivot_wider(names_from = overweight_obesity, values_from = dist) %>% 
    mutate(condition = 'overweight_obesity') %>% 
    bind_cols(p_val = round(t_test_macro("overweight_obesity", "energy")$p.val, 3))
```

    `summarise()` has grouped output by 'overweight_obesity'. You can override using the `.groups` argument.
    



```R
bind_rows(mets_macro, dm2_macro, overob_macro) 
```


<table class="dataframe">
<caption>A tibble: 9 × 5</caption>
<thead>
	<tr><th scope=col>source</th><th scope=col>FALSE</th><th scope=col>TRUE</th><th scope=col>condition</th><th scope=col>p_val</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>carbohydrate</td><td>54% 8%</td><td>55% 8%</td><td>met_s             </td><td>0.589</td></tr>
	<tr><td>lipids      </td><td>23% 5%</td><td>22% 5%</td><td>met_s             </td><td>0.841</td></tr>
	<tr><td>protein     </td><td>22% 5%</td><td>23% 5%</td><td>met_s             </td><td>0.394</td></tr>
	<tr><td>carbohydrate</td><td>55% 8%</td><td>55% 8%</td><td>dm2               </td><td>0.843</td></tr>
	<tr><td>lipids      </td><td>23% 5%</td><td>22% 6%</td><td>dm2               </td><td>0.214</td></tr>
	<tr><td>protein     </td><td>22% 5%</td><td>23% 5%</td><td>dm2               </td><td>0.537</td></tr>
	<tr><td>carbohydrate</td><td>55% 8%</td><td>55% 8%</td><td>overweight_obesity</td><td>0.954</td></tr>
	<tr><td>lipids      </td><td>23% 5%</td><td>23% 5%</td><td>overweight_obesity</td><td>0.450</td></tr>
	<tr><td>protein     </td><td>22% 5%</td><td>23% 5%</td><td>overweight_obesity</td><td>0.427</td></tr>
</tbody>
</table>




```R
bind_rows(energy_met_s, energy_dm2, energy_overob)
```


<table class="dataframe">
<caption>A tibble: 3 × 5</caption>
<thead>
	<tr><th scope=col>source</th><th scope=col>FALSE</th><th scope=col>TRUE</th><th scope=col>condition</th><th scope=col>p_val</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>energy</td><td>39442498</td><td>40132218</td><td>met_s             </td><td>0.675</td></tr>
	<tr><td>energy</td><td>39762398</td><td>39632204</td><td>d_diagnosed       </td><td>0.965</td></tr>
	<tr><td>energy</td><td>39372043</td><td>39852467</td><td>overweight_obesity</td><td>0.800</td></tr>
</tbody>
</table>




```R

```
