BEGIN TRANSACTION;

CREATE OR REPLACE TABLE metadata AS

--Filter the metadata
WITH filtered_metadata AS (
    SELECT *
    FROM metadata
    WHERE division IN ('Alaska', 'Washington', 'Oregon', 'Wyoming', 'Montana', 'Idaho', 'Massachusetts', 'Connecticut', 'Rhode Island', 'Vermont', 'New Hampshire', 'Maine')
      AND TRY_CAST(date AS DATE) IS NOT NULL
      AND TRY_CAST(date AS DATE) BETWEEN DATE '2020-01-01' AND DATE '2024-12-31'
),
--Count the size of each strata
group_counts AS (
    SELECT 
        division, 
        DATE_TRUNC('month', CAST(date AS DATE)) AS month,
        COUNT(*) AS group_size
    FROM filtered_metadata
    GROUP BY division, month
),
--Define the total size of the set
total_count AS (
    SELECT SUM(group_size) AS total_size FROM group_counts
),
--Allocate the number of sammples from each strata
allocations AS (
    SELECT 
        g.division,
        g.month,
        g.group_size,
        CEIL(g.group_size * 200000.0 / t.total_size) AS allocated_samples
    FROM group_counts g, total_count t
),
--Draw the rows
ranked_rows AS (
    SELECT 
        m.*,
        DATE_TRUNC('month', CAST(m.date AS DATE)) AS month,
        ROW_NUMBER() OVER (
            PARTITION BY m.division, DATE_TRUNC('month', CAST(m.date AS DATE))
            ORDER BY RANDOM()
        ) AS rn
    FROM filtered_metadata m
)
SELECT r.*
FROM ranked_rows r
JOIN allocations a 
  ON r.division = a.division AND r.month = a.month
WHERE r.rn <= a.allocated_samples;

--Subset pairs table from sampled sequences
CREATE OR REPLACE TABLE pairs AS
SELECT * FROM pairs
WHERE strain_1 IN (SELECT strain FROM metadata) 
AND strain_2 IN (SELECT strain FROM metadata);

COMMIT;
