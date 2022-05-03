1. Switch to your admin account #Skip if your user already has sudo permissions
   1. `su wag002adm` # replace with your admin account name
2. Start postgres
   1. `pg_ctl -D /opt/homebrew/var/postgres start`
   2. `ps aux | grep postgres` #to make sure that this process is running
3. Start psql and open database postgres, which is the database postgres uses itself to store roles, permissions, and structure:
   1. `psql postgres`
4. Create roles for the application, give login and CREATEDB permissions:
   1. `CREATE ROLE uta_admin WITH LOGIN CREATEDB;`
   2. `CREATE ROLE anonymous WITH LOGIN CREATEDB;`
5. List the users using the command 
   1. `\du`
6. Create the UTA Database object
   1. `CREATE DATABASE uta;`
7. Grant privileges to manage the database to uta_admin
   1. `GRANT ALL PRIVILEGES ON DATABASE uta TO uta_admin;`
8. Exit postgres
   1. `\q`
9. Download the UTA database and place it in the uta database object that you created before (**This step takes around 5 hours**). 
   1. `export UTA_VERSION=uta_20210129.pgd.gz\ncurl -O http://dl.biocommons.org/uta/$UTA_VERSION\ngzip -cdq ${UTA_VERSION} | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5432`
10. Set your UTA path
    1. `export UTA_DB_URL=postgresql://uta_admin@localhost:5432/uta/uta_20210129`
###Optional installation step
11. If you wanted to wait for the 5 hour update till later please follow these steps instead:
12. Download the UTA database and place it in the uta database object that you created before.
    1. `export UTA_VERSION=uta_20210129.pgd.gz
curl -O http://dl.biocommons.org/uta/$UTA_VERSION
gzip -cdq ${UTA_VERSION} | grep -v "^REFRESH MATERIALIZED VIEW" | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5432`
13. Run the refresh materialized view commands
    1. `REFRESH MATERIALIZED VIEW uta_20210129.exon_set_exons_fp_mv;`
    2. `REFRESH MATERIALIZED VIEW uta_20210129.tx_exon_set_summary_mv;`
    3. `REFRESH MATERIALIZED VIEW uta_20210129.tx_def_summary_mv;`
    4. `REFRESH MATERIALIZED VIEW uta_20210129.tx_similarity_mv;` #**This step will take 5 or more hours**