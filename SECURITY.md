# Security Policy

ProteinMutFlow is an alpha-stage local workflow orchestrator. It launches external scientific runtimes selected or discovered on the user's machine and therefore should be run only with trusted workflow files, structures, and mutation tables.

Do not include credentials, license files, proprietary structures, or private result data in a public issue. After the repository is published, use GitHub's private security-advisory channel for suspected vulnerabilities.

The project does not execute arbitrary shell fragments from `workflow.yaml`; backend paths and scientific selections are validated as data. This boundary should be preserved in contributions.
