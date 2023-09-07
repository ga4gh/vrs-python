1.  Switch to your admin account (Skip if your user already has `sudo` permissions):

        su wag002adm # replace with your admin account name`

2.  Start PostgreSQL:

    # OS and package manager dependent;

    # these commands and file locations assume an ARM Mac environment

    # and PostgreSQL installation with homebrew

    pg_ctl -D /opt/homebrew/var/postgres start
    ps aux | grep postgres # check that postgres processes are running

3.  Start `psql` and open the `postgres` database, which is the database PostgreSQL uses to store roles, permissions, and structure:

    psql postgres

4.  Create roles for the application, give login and CREATEDB permissions:

    CREATE ROLE uta_admin WITH LOGIN CREATEDB;
    CREATE ROLE anonymous WITH LOGIN CREATEDB;

5.  In `psql`, verify that the above users exist:

    \du

6.  In `psql`, create the UTA database:

    CREATE DATABASE uta;

7.  In `psql`, grant privileges to manage `uta` to `uta_admin`:

    GRANT ALL PRIVILEGES ON DATABASE uta TO uta_admin;

8.  Exit `psql`:

    \q

9.  Download the UTA database and place it in the `uta` PostgreSQL database that you created before (**This step takes around 5 hours** -- see optional instructions below for quick installation instructions):

    export UTA_VERSION=uta_20210129.pgd.gz
    curl -O http://dl.biocommons.org/uta/$UTA_VERSION
    gzip -cdq ${UTA_VERSION} | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5432

10. Set the UTA address environment variable (frequent users would be advised to add this to their login configs):

    export UTA_DB_URL=postgresql://uta_admin@localhost:5432/uta/uta_20210129

### Optional:

The majority of the load times described above come from refreshing materialized views. To delay this step, or avoid it entirely, you may optionally do the following instead:

If you wanted to wait for the 5 hour update till later please follow these steps instead:

9. Download the UTA database and place it in the `uta` PostgreSQL database that you created before:

   export UTA_VERSION=uta_20210129.pgd.gz
   curl -O http://dl.biocommons.org/uta/$UTA_VERSION
   gzip -cdq ${UTA_VERSION} | grep -v "^REFRESH MATERIALIZED VIEW" | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5432

10. Log back into `psql` and manually refresh the materialized views:

    REFRESH MATERIALIZED VIEW uta_20210129.exon_set_exons_fp_mv;
    REFRESH MATERIALIZED VIEW uta_20210129.tx_exon_set_summary_mv;
    REFRESH MATERIALIZED VIEW uta_20210129.tx_def_summary_mv;
    REFRESH MATERIALIZED VIEW uta_20210129.tx_similarity_mv; // This step will take 5 or more hours
