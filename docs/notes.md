Keys in record formats

`type`
  1. one-time-one-line
  2. one-time-multiple-lines
  3. multiple-time-one-line
  4. multiple-times-multiple-lines
  5. grouping
  6. other

`fields`: dictionary of base field names and type-byterange two-element lists

`continues`:  list of names of fields that is appended or extended by continuation records
 
`token_formats`:  dict keyed by field name listed in `fields`; keys:
  - `tokens`: dict keyed by token name; keys
    - `type`:  data type
    - `key`: optional string used to find token key in file if name is not used
  - `determinants`:  list of token names that group token:value pairs; a group must have same values of determinant tokens
  
`concatenate`:  dictionary that defines new fields formed by list-concatenation of base format fields; each key defines a new field whose value is formed by a list-like concatenation of base fields whose names are in the value

`mappings`: not used

`allowed`: dictionary of lists of allowed values for any given field

`determinants`: list of fields whose values determine whether this record is a new type or additional instance of a given type

`subrecords`: dictionary defining parsing of subrecords.  Keys:
  - `branchon`: Name of base field whose value determines which subrecord format is to be applied
  - `required`: (optional) boolean indicating whether we expect all base records to have a parseable branchon value
  - `formats`: dictionary of formats

`tables`: dictionary specifying how tables are parsed out of records