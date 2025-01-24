import flask
import numpy as np

blueprint = flask.Blueprint("ui", __name__)

@blueprint.route("/heartbeat", methods=["GET"])
def heartbeat():
    return "Up"

