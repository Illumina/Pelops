from ilmn.pelops.infrastructure import log_notifications


class TestLogNotificationFactory:
    def test_build(self):
        factory = log_notifications.LogNotificationServiceFactory()
        observed = factory.build()
        assert isinstance(observed, log_notifications.LogNotificationService)
