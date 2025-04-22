-- Use to trim the pair file down after cleaning the metadata
CREATE OR REPLACE TABLE pairs AS
SELECT * FROM pairs
    WHERE strain_1 IN (SELECT strain FROM metadata)
    AND strain_2 IN (SELECT strain FROM metadata);
