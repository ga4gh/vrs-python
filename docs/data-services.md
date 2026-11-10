# Set up SeqRepo and UTA

Some VRS operations need external data:

- [SeqRepo](https://github.com/biocommons/biocommons.seqrepo) provides reference
  sequences and aliases for validation, coordinate interpretation, and normalization.
- [UTA](https://github.com/biocommons/uta) provides transcript metadata and
  alignments for HGVS interpretation and mapping.

| Workflow | Data required |
| --- | --- |
| Model validation, serialization, identifier computation | Usually none; reference sequence validation requires SeqRepo |
| Allele normalization, SPDI/gnomAD translation, VCF annotation | Reference sequences and aliases, supplied by SeqRepo |
| HGVS transcript interpretation and mapping | Sequence data and, when transcript metadata or alignments are needed, UTA |

## Docker Compose

Docker Compose sets up SeqRepo REST and UTA, including the SeqRepo dataset.
You do not need a separate local SeqRepo installation for this setup.

Install `ga4gh.vrs[extras]` (see [installation](installation.md)) and
[Docker](https://docs.docker.com/get-started/get-docker/), then clone this repository:

```shell
git clone https://github.com/ga4gh/vrs-python.git
cd vrs-python
```

Run all Compose commands from the repository root, where `docker-compose.yml`
lives. Create the data volumes and start the services:

```shell
docker volume create uta_vol
docker volume create seqrepo_vol
docker compose up
```

Compose starts three containers:

| Service | Purpose |
| --- | --- |
| `seqrepo` | Downloads the `2024-12-20` dataset into `seqrepo_vol`, then exits |
| [seqrepo-rest-service](https://github.com/biocommons/seqrepo-rest-service) | Serves the downloaded sequences at `http://localhost:5000` |
| [uta](https://github.com/biocommons/uta) | Runs the `uta_20241220` transcript and alignment database on port `5432`, with data in `uta_vol` |

- The first startup may take 5–15 minutes, depending on your network and host.
- Allow enough disk space for the datasets. Later startups reuse the downloaded data.
- If SeqRepo REST starts before the download finishes, follow the
  [missing directory steps](#missing-seqrepo-directory).

These commands use Compose v2 (`docker compose`). If you still use v1, the command
is `docker-compose`; [upgrading to v2](https://docs.docker.com/compose/releases/migrate/)
is recommended.

### Configure Python clients

Set these environment variables so VRS-Python can connect to the services:

```shell
export GA4GH_VRS_DATAPROXY_URI='seqrepo+http://localhost:5000/seqrepo'
export UTA_DB_URL='postgresql://anonymous:anonymous@localhost:5432/uta/uta_20241220'
```

```python
from ga4gh.vrs.dataproxy import create_dataproxy

dp = create_dataproxy()  # Reads GA4GH_VRS_DATAPROXY_URI
print(dp.get_sequence("refseq:NM_000059.4", 0, 20))
```

The expected sequence is `AGAGGCGGAGCCGCTGTGGC`.

### Verify the services

In another terminal, check that `seqrepo-rest-service` and `uta` are running:

```shell
docker compose ps
```

Look for both services with status `Up` and ports `5000` and `5432` published.
The relevant columns should look similar to this (container names and elapsed
times may differ):

```text
NAME                                SERVICE                STATUS         PORTS
vrs-python-seqrepo-rest-service-1    seqrepo-rest-service   Up 2 minutes   0.0.0.0:5000->5000/tcp
vrs-python-uta-1                     uta                    Up 2 minutes   0.0.0.0:5432->5432/tcp
```

The `seqrepo` download container should exit after completing its work, so it
will not appear in the default `docker compose ps` output. Use
`docker compose ps --all` to see it; a successful download shows `Exited (0)`.

Test the two running services:

```shell
curl --write-out '\n' 'http://localhost:5000/seqrepo/1/sequence/refseq:NM_000059.4?end=20'
```

The sequence request should return `AGAGGCGGAGCCGCTGTGGC`.

Check UTA with the `psql` client:

```shell
psql -XAt 'postgresql://anonymous:anonymous@localhost:5432/uta' \
  -c 'select count(*) from uta_20241220.transcript'
```

The query should return a positive transcript count, such as `329090`.

## Local SeqRepo

A local sequence store avoids REST requests. Use this as an alternative to the
Compose-managed SeqRepo volume. Install SeqRepo and prepare its default data
directory:

```shell
pip install biocommons.seqrepo
sudo mkdir -p /usr/local/share/seqrepo
sudo chown "$USER" /usr/local/share/seqrepo
```

Choose a dataset listed by `seqrepo list-remote-instances`, then download it:

```shell
seqrepo list-remote-instances
export SEQREPO_VERSION=2024-12-20  # Or a newer available dataset
seqrepo pull -i "$SEQREPO_VERSION"
seqrepo update-latest
export GA4GH_VRS_DATAPROXY_URI='seqrepo+file:///usr/local/share/seqrepo/latest'
```

See [SeqRepo's installation documentation](https://github.com/biocommons/biocommons.seqrepo)
for alternate storage locations.

If a completed download raises a `PermissionError` when renaming its temporary
directory, move it manually. For example:

```shell
sudo mv /usr/local/share/seqrepo/2024-12-20._fkuefgd /usr/local/share/seqrepo/2024-12-20
seqrepo update-latest
```

Replace those paths with the temporary and final paths from your error message.

## Local UTA

To install UTA locally, follow the instructions in the
[UTA repository](https://github.com/biocommons/uta).

## Configuration

| Setting | Purpose |
| --- | --- |
| `GA4GH_VRS_DATAPROXY_URI` | Read automatically by `create_dataproxy()` when no explicit `uri` is supplied |
| `UTA_DB_URL` | HGVS library's UTA connection, including the schema name |
| `SEQREPO_URI` | Read by some notebooks and examples and passed as an explicit `uri`, overriding `GA4GH_VRS_DATAPROXY_URI`; uses the same `seqrepo+` format |

Data-proxy URIs must include the `seqrepo+` provider prefix. Supported formats:

- Local absolute path: `seqrepo+file:///usr/local/share/seqrepo/latest`
- Local relative path: `seqrepo+:../relative/path/to/seqrepo/latest`
- HTTP service: `seqrepo+http://example.org:5000/seqrepo`
- HTTPS service: `seqrepo+https://example.org/seqrepo`

You must supply a URI or set `GA4GH_VRS_DATAPROXY_URI`.

## Troubleshooting

### Missing SeqRepo directory

If `seqrepo-rest-service` exits with `OSError: Unable to open SeqRepo directory
/usr/local/share/seqrepo/2024-12-20`, check the download container's logs:

```shell
docker compose logs seqrepo
```

Wait for the download to finish, then start the REST service again:

```shell
docker compose up seqrepo-rest-service
```

You can start services individually this way to isolate a startup failure. To
inspect the REST service or UTA logs, run:

```shell
docker compose logs seqrepo-rest-service uta
```

### Port already in use

Another PostgreSQL server may already use port `5432`. A UTA connection error
such as `FATAL: role "anonymous" does not exist` can mean your client reached that
server. Another process may also occupy SeqRepo REST's port `5000`.

Check the port with `sudo lsof -i :5432` or `sudo lsof -i :5000`. Choose an unused
port; `sudo lsof -i :5433`, for example, should produce no output if it is free.

In `docker-compose.yml`, change the first number in the service's port mapping.
The first number is the host port; the second is the container port. For example,
use `5433:5432` for UTA or `5001:5000` for SeqRepo REST:

```yaml
# Under the uta service:
ports:
  - 5433:5432

# Under the seqrepo-rest-service service:
ports:
  - 5001:5000
```

Run `docker compose up` again and check the new ports with `lsof`. Update the
environment variables to use the new ports:

```shell
export UTA_DB_URL='postgresql://anonymous:anonymous@localhost:5433/uta/uta_20241220'
export GA4GH_VRS_DATAPROXY_URI='seqrepo+http://localhost:5001/seqrepo'
export SEQREPO_URI="$GA4GH_VRS_DATAPROXY_URI"
```

Test using the new host ports:

```shell
psql -XAt 'postgresql://anonymous:anonymous@localhost:5433/uta' \
  -c 'select count(*) from uta_20241220.transcript'
curl --write-out '\n' 'http://localhost:5001/seqrepo/1/sequence/refseq:NM_000059.4?end=20'
```

### Other connection and data errors

- **UTA role or schema missing:** check that the URL reaches the intended UTA
  container, rather than another PostgreSQL server, and names the installed schema.
- **Unknown sequence or alias:** check accession versions, assembly, and whether
  the dataset contains the requested reference.
- **Invalid data-proxy scheme:** include `seqrepo+`; a bare HTTP URL is insufficient.
