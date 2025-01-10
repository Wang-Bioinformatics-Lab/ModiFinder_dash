build:
	docker build -t gnpslibrary . 

bash:
	docker run -it --rm gnpslibrary /bin/bash

server-compose-build-nocache:
	docker-compose --compatibility build --no-cache

server-compose-interactive:
	docker-compose --compatibility build
	docker-compose --compatibility up

server-compose:
	docker-compose --compatibility build
	docker-compose --compatibility up -d

server-compose-production:
	docker-compose --compatibility build
	docker-compose --compatibility -f docker-compose.yml -f docker-compose-prod.yml up -d

server-compose-production-interative:
	docker-compose --compatibility build
	docker-compose --compatibility -f docker-compose.yml -f docker-compose-prod.yml up

attach:
	docker exec -i -t mod-site /bin/bash
