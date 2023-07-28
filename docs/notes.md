Keys in record formats

- type

  Integer 1, 2, 3, 4, 5, or 6
  - 1: one-time-one-line
  - 2: one-time-multiple-lines
  - 3: multiple-time-one-line
  - 4: multiple-times-multiple-lines
  - 5: grouping
  - 6: other
- fields

  dictionary of base field names and type-byterange two-element lists

- continues

  list of names of fields that is appended or extended by continuation records
 
- tokens

  dictionary keyed by field name; values are dictionary of token-dicts keyed on token names.  A token-dict has a type element and an optional `associated_to` element for grouping

- concatenate

  dictionary that defines new fields formed by list-concatenation of base format fields

- mappings

    not used

- allowed

    dictionary of lists of allowed values for any given field

- determinants

    list of fields whose values determine whether this record is a new type or additional instance of a given type

- subrecords

    the record line can be declared as a subrecord via the `branchon` key whose value is the name of the field that determines the subrecord format (dictionary `formats`) used to parse the record line.